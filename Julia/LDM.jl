# ==============================
# Formulates the problem and solves it
# ==============================
using JuMP
using PyPlot
using Gurobi, KNITRO
include("utilities.jl")

# ==============================
# Formulates the problem and solves it
# ==============================
function GenerateInputs(types, fm)
    ntypes = length(types)
    nsupp = length(fm)
    Theta = combwithrep(nsupp, ntypes)

    # joint distribution based on marginals
    f = Dict()
    for perm in Theta
        f[perm] = prod([fm[i][perm[i]] for i=1:nsupp])
    end

    # ------------------------------------------------------------------------------
    # CENTRALIZED WITH IC AND IR
    # ------------------------------------------------------------------------------
    nvars = length(Theta)*nsupp # here we include transfers and allocations
    sts = length(Theta) # size of type space

    # this dict tell us the probability of all other subjects being of their type $f_{-i}(\theta_{-i})$
    f_woi = Dict(i=>Dict(j=>0.0 for j in Theta) for i in 1:nsupp)

    for i in keys(f_woi)
        supp_woi = [j for j in 1:nsupp if j!=i]
        for perm in Theta
            f_woi[i][perm] = prod([fm[j][perm[j]] for j in supp_woi])
        end
    end

    if nsupp == 1
        f_woi = Dict(i=>Dict(j=>1.0 for j in Theta) for i in 1:nsupp)
    end

    # individual rationality constraints
    bIR_x = zeros((ntypes*nsupp, nvars))
    bIR_t = zeros((ntypes*nsupp, nvars))
    row = 1
    for i in 1:nsupp
        for j in 1:ntypes
            # we assume that the type of employee i is j
            # now we find all scenarios where employee i is of type j
            idxs = [k for k in 1:length(Theta) if Theta[k][i] == j]
            for k in idxs
                bIR_x[row, nsupp*(k-1)+i] = -types[j]*f_woi[i][Theta[k]]
                bIR_t[row, nsupp*(k-1)+i] = f_woi[i][Theta[k]]
            end
            row=row+1
        end
    end

    bIC_x = zeros((ntypes*(ntypes-1)*nsupp,nvars))
    bIC_t = zeros((ntypes*(ntypes-1)*nsupp,nvars))
    row = 1
    for i in 1:nsupp
        for j in 1:ntypes
            other_types = [k for k in 1:ntypes if k!=j]
            for k in other_types
                # print j, other_types, row
                idxs_j = [l for l in 1:length(Theta) if Theta[l][i] == j]
                idxs_k = [l for l in 1:length(Theta) if Theta[l][i] == k]
                for l in idxs_j
                    bIC_x[row, nsupp*(l-1)+i] = -types[j]*f_woi[i][Theta[l]]
                    bIC_t[row, nsupp*(l-1)+i] = f_woi[i][Theta[l]]
                end
                for l in idxs_k
                    bIC_x[row, nsupp*(l-1)+i] = types[j]*f_woi[i][Theta[l]]
                    bIC_t[row, nsupp*(l-1)+i] = -f_woi[i][Theta[l]]
                end
                row+=1
            end
        end
    end

    bG_x = vcat(bIR_x, bIC_x)
    bG_t = vcat(bIR_t, bIC_t)

    bA = zeros((length(Theta), nvars))
    bh = zeros(size(bG_x)[1])
    bb = ones(length(Theta))

    wq_t = zeros(nvars)
    q_t = zeros(nvars)

    for i in 1:length(Theta)
        bA[i,nsupp*(i-1)+1:nsupp*i] = ones(nsupp)
        wq_t[nsupp*(i-1)+1:nsupp*i] = f[Theta[i]]*ones(nsupp)
        q_t[nsupp*(i-1)+1:nsupp*i] = ones(nsupp)
    end

    return nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t
end

function InputsObjectiveLM(nalpha, nGamma, nsupp, nvars, sts)
    nD = inv(nGamma)
    nc = nD*nalpha
    D = zeros((nvars, nvars)) # matrix D in diagonal
    wD = zeros((nvars, nvars)) # matrix D weighted by distribution of types
    wq_x = zeros(nvars)
    q_x = zeros(nvars)
    for i in 1:sts
        D[nsupp*(i-1)+1:nsupp*i,nsupp*(i-1)+1:nsupp*i] = nD
        wD[nsupp*(i-1)+1:nsupp*i,nsupp*(i-1)+1:nsupp*i] = f[Theta[i]]*nD
        wq_x[nsupp*(i-1)+1:nsupp*i] = f[Theta[i]]*nc
        q_x[nsupp*(i-1)+1:nsupp*i] = nc
    end
    return D, wD, q_x, wq_x
end

function CheckFeasibility(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, nalpha, nGamma, version, elastic=false, x0=nothing, t0=nothing)
    D, wD, c_x, wc_x = InputsObjectiveLM(nalpha, nGamma, nsupp, nvars, sts)

    if version == "centralized"
        m = Model(solver=GurobiSolver(Presolve=0, MIPGap = 1e-12))
    elseif version == "decentralized"
        m = Model(solver=KnitroSolver(mip_method = KTR_MIP_METHOD_BB, honorbnds=0,
                                      ms_enable = 1, ms_maxsolves = 5000,
                                      algorithm = KTR_ALG_ACT_CG,
                                      outmode = KTR_OUTMODE_SCREEN,
                                      KTR_PARAM_OUTLEV = KTR_OUTLEV_ALL,
                                      KTR_PARAM_MIP_OUTINTERVAL = 1,
                                      KTR_PARAM_MIP_MAXNODES = 10000,
                                      KTR_PARAM_HESSIAN_NO_F = KTR_HESSIAN_NO_F_ALLOW))
    else
        println("***ERROR: unknown version")
        exit()
    end

    if x0 != nothing
        @variable(m, x[i=1:nvars]==x0[i])
        @variable(m, t[i=1:nvars]==t0[i])
    else
        @variable(m, x[1:nvars]>=0)
        @variable(m, t[1:nvars]>=0)
    end
    @variable(m, z)

    # constraints centralized version: IR, IC, feas
    @constraint(m, bG_x*x + bG_t*t .>= bh) # IR + IC
    if elastic == false
        @constraint(m, bA*x .== bb) # feas
    end
    @constraint(m, z <=  wc_x'*x - wq_t'*t - 0.5*x'*wD*x )

    if version == "decentralized"
        @variable(m, p[1:nvars]>=0)
        @variable(m, u[1:nvars]>=0)
        @variable(m, v[1:nvars])
        @constraint(m, kkt_opt, p - c_x + D*x - u + v .== 0)
        @NLconstraint(m, var_def[i in 1:nvars], t[i]-x[i]*p[i] == 0)
        @NLconstraint(m, kkt_comp[i in 1:nvars], x[i]*u[i] == 0)
        @constraint(m, kkt_cons[i in 1:sts], v[nsupp*(i-1)+1]-v[nsupp*(i-1)+2] == 0)
    end
    @objective(m, Max, 1)
    print(m)
    status = solve(m)
    obj = -getobjectivevalue(m)
    x_vals = getvalue(x)
    t_vals = getvalue(t)
    transfers = wq_t'*t_vals
    println("===============================")
    println("Objective value: ", getobjectivevalue(m))
    println("Allocations: ", getvalue(x))
    println("Transfers: ", getvalue(t))
    println("===============================")
    return obj, transfers, x_vals, t_vals
end

function SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, nalpha, nGamma, version, elastic=false, x0=nothing, t0=nothing)
    D, wD, c_x, wc_x = InputsObjectiveLM(nalpha, nGamma, nsupp, nvars, sts)

    if version == "centralized"
        m = Model(solver=GurobiSolver(Presolve=0, MIPGap = 1e-12))
    elseif version == "decentralized"
        m = Model(solver=KnitroSolver(mip_method = KTR_MIP_METHOD_BB, honorbnds=0,
                                      ms_enable = 1, ms_maxsolves = 5000,
                                      algorithm = KTR_ALG_ACT_CG,
                                      outmode = KTR_OUTMODE_SCREEN,
                                      KTR_PARAM_OUTLEV = KTR_OUTLEV_ALL,
                                      KTR_PARAM_MIP_OUTINTERVAL = 1,
                                      KTR_PARAM_MIP_MAXNODES = 10000,
                                      KTR_PARAM_HESSIAN_NO_F = KTR_HESSIAN_NO_F_ALLOW))
    else
        println("***ERROR: unknown version")
        exit()
    end

    if x0 != nothing
        @variable(m, x[i=1:nvars]==x0[i])
        @variable(m, t[i=1:nvars]==t0[i])
    else
        @variable(m, x[1:nvars]>=0)
        @variable(m, t[1:nvars]>=0)
    end
    @variable(m, z)

    # constraints centralized version: IR, IC, feas
    @constraint(m, bG_x*x + bG_t*t .>= bh) # IR + IC
    if elastic == false
        @constraint(m, bA*x .== bb) # feas
    end
    @constraint(m, z <=  wc_x'*x - wq_t'*t - 0.5*x'*wD*x )

    if version == "decentralized"
        @variable(m, p[1:nvars]>=0)
        @variable(m, u[1:nvars]>=0)
        @variable(m, v[1:nvars])
        @constraint(m, kkt_opt, p - c_x + D*x - u + v .== 0)
        @NLconstraint(m, var_def[i in 1:nvars], t[i]-x[i]*p[i] == 0)
        @NLconstraint(m, kkt_comp[i in 1:nvars], x[i]*u[i] == 0)
        @constraint(m, kkt_cons[i in 1:sts], v[nsupp*(i-1)+1]-v[nsupp*(i-1)+2] == 0)
    end

    @objective(m, Max, z)
    print(m)
    status = solve(m)
    obj = -getobjectivevalue(m)
    x_vals = getvalue(x)
    t_vals = getvalue(t)
    transfers = wq_t'*t_vals
    println("===============================")
    println("Objective value: ", getobjectivevalue(m))
    println("Allocations: ", getvalue(x))
    println("Transfers: ", getvalue(t))
    println("===============================")
    return obj, transfers, x_vals, t_vals
end


function SimulateOptimization(types, fm, loc, gammas, version, elastic=false)
    println("Generating Inputs")
    nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t = GenerateInputs(types, fm)
    objs, trans = [],[]
    for i=1:length(gammas)
        nGamma = gammas[i]
        println("Solving the problem for ", delta)
        if version == "decentralized"
            # obtain initial solution
            obj_cent, tra_cent, x_cent, t_cent = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, nalpha, nGamma, "centralized", elastic)
            obj, tra, xs, ts = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, nalpha, nGamma, version, elastic, x_cent, t_cent)
        else
            obj, tra, xs, ts = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, nalpha, nGamma, version, elastic)
        end
        push!(objs, obj)
        push!(trans, tra)
    end
    println(objs)
    println(trans)
    plot(deltas, objs, color="blue", linewidth=2.0, linestyle="--")
    plot(deltas, trans, color="red", linewidth=2.0, linestyle="--")
    if elastic
        outfile = string("plots/objective_and_transfers_LDM_elastic_", version, ".pdf")
    else
        outfile = string("plots/objective_and_transfers_LDM_inelastic_", version, ".pdf")
    end

    savefig(outfile)
    show()
end

function FormulateAndSolve(types, fm, nalpha, nGamma, version, elastic=false)

    ntypes = length(types)
    nsupp = length(fm)
    Theta = combwithrep(nsupp, ntypes)

    # joint distribution based on marginals
    f = Dict()
    for perm in Theta
        f[perm] = prod([fm[i][perm[i]] for i=1:nsupp])
    end



    # bulding cumulative distribution

    # ------------------------------------------------------------------------------
    # CENTRALIZED WITH IC AND IR
    # ------------------------------------------------------------------------------
    nvars = length(Theta)*nsupp # here we include transfers and allocations
    sts = length(Theta) # size of type space



    # this dict tell us the probability of all other subjects being of their type $f_{-i}(\theta_{-i})$
    f_woi = Dict(i=>Dict(j=>0.0 for j in Theta) for i in 1:nsupp)
    # println(f_woi)
    # println(Theta)


    for i in keys(f_woi)
        supp_woi = [j for j in 1:nsupp if j!=i]
        for perm in Theta
            f_woi[i][perm] = prod([fm[j][perm[j]] for j in supp_woi])
        end
    end

    if nsupp == 1
        f_woi = Dict(i=>Dict(j=>1.0 for j in Theta) for i in 1:nsupp)
    end

    # individual rationality constraints
    bIR_x = zeros((ntypes*nsupp, nvars))
    bIR_t = zeros((ntypes*nsupp, nvars))
    row = 1
    for i in 1:nsupp
        for j in 1:ntypes
            # we assume that the type of employee i is j
            # now we find all scenarios where employee i is of type j
            idxs = [k for k in 1:length(Theta) if Theta[k][i] == j]
            for k in idxs
                bIR_x[row, nsupp*(k-1)+i] = -types[j]*f_woi[i][Theta[k]]
                bIR_t[row, nsupp*(k-1)+i] = f_woi[i][Theta[k]]
            end
            row=row+1
        end
    end
    # println(bIR_x)
    # println(bIR_t)


    # # checking IR constraints
    # println(fm)
    # println(types)
    # println(f_woi)
    # for i in 1:size(bIR_x)[1]
    #     println(bIR_x[i,:], " ", bIR_t[i,:])
    # end
    # exit()



    bIC_x = zeros((ntypes*(ntypes-1)*nsupp,nvars))
    bIC_t = zeros((ntypes*(ntypes-1)*nsupp,nvars))
    row = 1
    for i in 1:nsupp
        for j in 1:ntypes
            other_types = [k for k in 1:ntypes if k!=j]
            for k in other_types
                # print j, other_types, row
                idxs_j = [l for l in 1:length(Theta) if Theta[l][i] == j]
                idxs_k = [l for l in 1:length(Theta) if Theta[l][i] == k]
                for l in idxs_j
                    bIC_x[row, nsupp*(l-1)+i] = -types[j]*f_woi[i][Theta[l]]
                    bIC_t[row, nsupp*(l-1)+i] = f_woi[i][Theta[l]]
                end
                for l in idxs_k
                    bIC_x[row, nsupp*(l-1)+i] = types[j]*f_woi[i][Theta[l]]
                    bIC_t[row, nsupp*(l-1)+i] = -f_woi[i][Theta[l]]
                end
                row+=1
            end
        end
    end

    # checking IC constraints
    # println(fm)
    # println(types)
    # println(f_woi)
    # for i in 1:size(bIC_x)[1]
    #     println(bIR_x[i,:], " ", bIR_t[i,:])
    #     println(bIC_x[i,:], " ", bIC_t[i,:])
    #     println("-------------------------------------------------")
    # end
    # exit()


    bG_x = vcat(bIR_x, bIC_x)
    bG_t = vcat(bIR_t, bIC_t)
    # println(bG)
    # println(size(bG_x))
    # println(size(bG_t))

    nD = inv(nGamma)
    nc = nD*nalpha
    D = zeros((nvars, nvars)) # matrix D in diagonal
    wD = zeros((nvars, nvars)) # matrix D weighted by distribution of types

    bA = zeros((length(Theta), nvars))
    bh = zeros(size(bG_x)[1])
    bb = ones(length(Theta))

    wq_x = zeros(nvars)
    wq_t = zeros(nvars)
    q_x = zeros(nvars)
    q_t = zeros(nvars)


    for i in 1:length(Theta)
        D[nsupp*(i-1)+1:nsupp*i,nsupp*(i-1)+1:nsupp*i] = nD
        wD[nsupp*(i-1)+1:nsupp*i,nsupp*(i-1)+1:nsupp*i] = f[Theta[i]]*nD
        bA[i,nsupp*(i-1)+1:nsupp*i] = ones(nsupp)
        wq_x[nsupp*(i-1)+1:nsupp*i] = f[Theta[i]]*nc
        wq_t[nsupp*(i-1)+1:nsupp*i] = f[Theta[i]]*ones(nsupp)
        q_x[nsupp*(i-1)+1:nsupp*i] = nc
        q_t[nsupp*(i-1)+1:nsupp*i] = ones(nsupp)
    end

    # checking feasibility constraints
    # for i in 1:size(bA)[1]
    #     println(bA[i,:], " ", bb[i])
    # end

    # checking matrix D
    # for i in 1:size(D)[1]
    #     println(D[i,:])
    # end
    # println(f)
    # for i in 1:size(wD)[1]
    #     println(wD[i,:])
    # end

    # println(wq_x)
    # println(q_x)
    # println(wq_t)
    # println(q_t)
    #
    # exit()
    # x0 = [0.363636, 0.636364, 0.672727, 0.327273, 0.634091, 0.365909, 0.131818, 0.868182, 0.440909, 0.559091, 0.402273, 0.597727, 0.131818, 0.868182, 0.440909, 0.559091, 0.402273, 0.597727]
    # t0 = [0.1,0.1,0.1,0.5,0.1,0.45,0.4,0.1,0.4,0.5,0.4,0.45,0.4,0.1,0.4,0.5,0.4,0.45]

    # f_0 = wq_x'*x0 - wq_t'*t0 - 0.5*x0'*wD*x0
    # println(f_0)
    # println(bA*x0)
    # println(bG_x*x0 + bG_t*t0 )
    # # 0.004353273328039364
    # exit()

    # OPTIMIZATION
    if version == "centralized"
        m = Model(solver=GurobiSolver(Presolve=0, MIPGap = 1e-12))
    elseif version == "decentralized"
        m = Model(solver=KnitroSolver(mip_method = KTR_MIP_METHOD_BB,
                                      ms_enable = 1, ms_maxsolves = 1000,
                                      algorithm = KTR_ALG_ACT_CG,
                                      outmode = KTR_OUTMODE_SCREEN,
                                      KTR_PARAM_OUTLEV = KTR_OUTLEV_ALL,
                                      KTR_PARAM_MIP_OUTINTERVAL = 1,
                                      KTR_PARAM_MIP_MAXNODES = 10000,
                                      KTR_PARAM_HESSIAN_NO_F = KTR_HESSIAN_NO_F_ALLOW))
    else
        println("***ERROR: unknown version")
        exit()
    end

    @variable(m, x[1:nvars]>=0)
    @variable(m, t[1:nvars]>=0)
    @variable(m, z)

    # constraints centralized version: IR, IC, feas
    @constraint(m, bG_x*x + bG_t*t .>= bh) # IR + IC
    @constraint(m, bA*x .== bb) # feas
    @constraint(m, z <=  wq_x'*x - wq_t'*t - 0.5*x'*wD*x )

    if version == "decentralized"
        @variable(m, p[1:nvars]>=0)
        @variable(m, u[1:nvars]>=0)
        @variable(m, v[1:nvars])
        @constraint(m, kkt_opt, p - q_x + D*x - u + v .== 0)
        @NLconstraint(m, var_def[i in 1:nvars], t[i]-x[i]*p[i] == 0)
        @NLconstraint(m, kkt_comp[i in 1:nvars], x[i]*u[i] == 0)
        @constraint(m, kkt_cons[i in 1:sts], v[nsupp*(i-1)+1]-v[nsupp*(i-1)+2] == 0)
    end

    @objective(m, Max, z)

    print(m)

    status = solve(m)

    println("Objective value: ", getobjectivevalue(m))

    println("Allocations: ", getvalue(x))
end


# IMPORTANT: types must be sorted in increasing order
### INPUTS ###
# types = Dict(1=>0.1, 2=>0.2)
# fm = Dict(1=>[0.6,0.4],2=>[0.75,0.25])
types = Dict(1=>0.1, 2=>0.2, 3=>0.25)
fm = Dict(1=>[0.25,0.5, 0.25],2=>[0.2, 0.6, 0.2])
V = check_vc_increasing(fm)

# alphas
a_1, a_2 = 0.25, 0.5
# matrix for demand
r_1, r_2 = 1,1.1
gamma = 0.5

if r_1 + r_2 < 2*gamma
    println("***ERROR: does not satisfy consistency check")
    exit()
end

# build input matrices
if length(fm) == 2
    nalpha = [a_1; a_2]
    nGamma = [r_1 -gamma;-gamma r_2]
end
if length(fm) == 1
    nalpha = [a_1]
    nGamma = [r_1]
end

version = "decentralized" # or decentralized

FormulateAndSolve(types, fm, nalpha, nGamma, version)

exit()
