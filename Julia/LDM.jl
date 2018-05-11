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

    return nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, f, Theta
end

function GenerateExPostInputs(types, Theta, nsupp, nvars)
    ntypes = length(types)
    # individual rationality constraints
    sts = length(Theta)
    bIR_x = zeros((sts*nsupp, nvars))
    bIR_t = zeros((sts*nsupp, nvars))
    row = 1

    for i in 1:nsupp
        for j in 1:ntypes
            # we assume that the type of employee i is j
            # now we find all scenarios where employee i is of type j
            idxs = [k for k in 1:length(Theta) if Theta[k][i] == j]
            for k in idxs
                bIR_x[row, nsupp*(k-1)+i] = -types[j]
                bIR_t[row, nsupp*(k-1)+i] = 1
                row=row+1
            end
        end
    end
    bh = zeros(size(bIR_x)[1])

    return bIR_x, bIR_t, bh
end

function InputsObjectiveLM(nalpha, nGamma, f, Theta, nsupp, ntypes, nvars, sts)
    nD = inv(nGamma)
    nc = nD*nalpha

    D = zeros((nvars, nvars)) # matrix D in diagonal
    wD = zeros((nvars, nvars)) # matrix D weighted by distribution of types
    wq_x = zeros(nvars)
    q_x = zeros(nvars)
    for i in 1:sts
        # unweighted parameters for objective function
        D[nsupp*(i-1)+1:nsupp*i,nsupp*(i-1)+1:nsupp*i] = nD
        q_x[nsupp*(i-1)+1:nsupp*i] = nc
        # weighted by probabilities of each scenario
        wD[nsupp*(i-1)+1:nsupp*i,nsupp*(i-1)+1:nsupp*i] = f[Theta[i]]*nD
        wq_x[nsupp*(i-1)+1:nsupp*i] = f[Theta[i]]*nc
    end
    return D, wD, q_x, wq_x
end

function CheckFeasibility(nsupp, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, D, wD, c_x, wc_x, epG_x, epG_t, eph, version, elastic=false, expostir=false, x0=nothing, t0=nothing)

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

function SolveOptimization(nsupp, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, D, wD, c_x, wc_x, epG_x, epG_t, eph, version, elastic=false, expostir=false, x0=nothing, t0=nothing)

    if version == "centralized"
        m = Model(solver=GurobiSolver(Presolve=0, MIPGap = 1e-12))
    elseif version == "decentralized"
        m = Model(solver=KnitroSolver(mip_method = KTR_MIP_METHOD_BB, honorbnds=0,
                                      ms_enable = 1, ms_maxsolves = 500,
                                      algorithm = KTR_ALG_ACT_CG,
                                      outmode = KTR_OUTMODE_SCREEN,
                                      KTR_PARAM_OUTLEV = KTR_OUTLEV_ALL,
                                      KTR_PARAM_MIP_OUTINTERVAL = 1,
                                      KTR_PARAM_MIP_MAXNODES = 10000,
                                      KTR_PARAM_HESSIAN_NO_F = KTR_HESSIAN_NO_F_ALLOW)) #par_numthreads = 2, par_msnumthreads = 1,
    else
        println("***ERROR: unknown version")
        exit()
    end

    if x0 != nothing
        @variable(m, x[i=1:nvars]>=0, start=x0[i])
        @variable(m, t[i=1:nvars]>=0, start=t0[i])
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

    # special constraints
    if expostir == true
        @constraint(m, epG_x*x + epG_t*t .>= eph) # IR + IC
    end
    if elastic == false
        @constraint(m, bA*x .== bb) # feas
    end

    if version == "decentralized"
        @variable(m, p[1:nvars]>=0)
        @variable(m, u[1:nvars]>=0)
        if elastic == false
            @variable(m, v[1:nvars])
            @constraint(m, kkt_cons[i in 1:sts], v[nsupp*(i-1)+1]-v[nsupp*(i-1)+2] == 0)
            @constraint(m, kkt_opt, p - c_x + D*x - u + v .== 0)
        else
            @constraint(m, kkt_opt, p - c_x + D*x - u .== 0)
        end
        @NLconstraint(m, var_def[i in 1:nvars], t[i]-x[i]*p[i] == 0)
        @NLconstraint(m, kkt_comp[i in 1:nvars], x[i]*u[i] == 0)
    end

    @objective(m, Max, z)
    print(m)
    status = solve(m)
    obj = getobjectivevalue(m)
    x_vals = getvalue(x)
    t_vals = getvalue(t)
    println(c_x)
    transfers = wq_t'*t_vals
    println("===============================")
    println("Objective value: ", getobjectivevalue(m))
    println("Allocations: ", getvalue(x))
    println("Transfers: ", getvalue(t))
    println("===============================")
    return obj, transfers, x_vals, t_vals
end

function FormulateAndSolve(types, fm, nalpha, nGamma, elastic=false, expostir=false)
    nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, distr, Theta = GenerateInputs(types, fm)
    D, wD, c_x, wc_x = InputsObjectiveLM(nalpha, nGamma, distr, Theta, nsupp, ntypes, nvars, sts)
    epG_x, epG_t, eph = GenerateExPostInputs(types, Theta, nsupp, nvars)
    obj_cent, transfers_cent, x_cent, t_cent = SolveOptimization(nsupp, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, D, wD, c_x, wc_x, epG_x, epG_t, eph, "centralized", elastic, expostir)
    obj_dec, transfers_dec, x_dec, t_dec = SolveOptimization(nsupp, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, D, wD, c_x, wc_x, epG_x, epG_t, eph, "decentralized", elastic, expostir, x_cent, t_cent )
    return obj_cent, transfers_cent, x_cent, t_cent, obj_dec, transfers_dec, x_dec, t_dec
end
