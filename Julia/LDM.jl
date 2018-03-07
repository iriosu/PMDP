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

    return nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, f, Theta
end

function InputsObjectiveLM(nalpha, nGamma, f, Theta, nsupp, ntypes, nvars, sts)
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

function CheckFeasibility(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, nalpha, nGamma, distr, Theta, version, elastic=false, x0=nothing, t0=nothing)
    D, wD, c_x, wc_x = InputsObjectiveLM(nalpha, nGamma, distr, Theta, nsupp, ntypes, nvars, sts)

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

function SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, nalpha, nGamma, distr, Theta, version, elastic=false, x0=nothing, t0=nothing)
    D, wD, c_x, wc_x = InputsObjectiveLM(nalpha, nGamma, distr, Theta, nsupp, ntypes, nvars, sts)

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
                                      KTR_PARAM_HESSIAN_NO_F = KTR_HESSIAN_NO_F_ALLOW)) #par_numthreads = 2, par_msnumthreads = 1,
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
    println(c_x)
    transfers = wq_t'*t_vals
    println("===============================")
    println("Objective value: ", getobjectivevalue(m))
    println("Allocations: ", getvalue(x))
    println("Transfers: ", getvalue(t))
    println("===============================")
    return obj, transfers, x_vals, t_vals
end

function SimulateOptimization(types, fm, a_j, gamma_ii, gamma_ij, elastic=false)
    nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, distr, Theta = GenerateInputs(types, fm)
    if elastic
        f = open("simulations_outcome_elastic.txt", "w")
    else
        f = open("simulations_outcome_inelastic.txt", "w")
    end
    for i1=1:length(a_j)
        # to start assume that qualities are the same for both suppliers and same price sensitivities
        for i2=i1:length(a_j)
            alpha = [a_j[i1]; a_j[i2]]
            for j1=1:length(gamma_ii)
                for j2=j1:length(gamma_ii)
                    for k=1:length(gamma_ij)
                        Gamma = [gamma_ii[j1] -gamma_ij[k]; -gamma_ij[k] gamma_ii[j2]]
                        if prod([check_diagonal_dominant(Gamma), check_consistency_demand_matrix(Gamma), isposdef(Gamma)]) == false
                            println("***WARNING: the matrix does not satisfy the assumptions")
                            continue
                        end
                        obj_cent, transfers_cent, x_cent, t_cent = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, alpha, Gamma, distr, Theta, "centralized", elastic)
                        obj_dec, transfers_dec, x_dec, t_dec = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, alpha, Gamma, distr, Theta, "decentralized", elastic, x_cent, t_cent)
                        # check if set of active suppliers is the same
                        boo = true
                        for l=1:length(x_cent)
                            if x_cent[l] > 1e-5 && x_dec[l] < 1e-5
                                boo = false
                                break
                            elseif x_cent[l] < 1e-5 && x_dec[l] > 1e-5
                                boo = false
                                break
                            else
                                continue
                            end
                        end

                        outstr = string(a_j[i1], ";", a_j[i2], ";", gamma_ii[j1], ";", gamma_ii[j2], ";", gamma_ij[k], ";",
                                        obj_cent, ";", obj_dec, ";", boo)
                        write(f, "$outstr \n")
                        exit()
                    end
                end
            end
        end
    end
    close(f)
end

function FormulateAndSolve(types, fm, nalpha, nGamma, version, elastic=false)
    nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, distr, Theta = GenerateInputs(types, fm)
    obj, transfers, x_vals, t_vals = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, nalpha, nGamma, distr, Theta, version, elastic)
end



# IMPORTANT: types must be sorted in increasing order
### INPUTS ###
types = Dict(1=>10, 2=>12)
fm = Dict(1=>[0.6,0.4],2=>[0.5,0.5])
# types = Dict(1=>0.1, 2=>0.2, 3=>0.25)
# fm = Dict(1=>[0.25,0.5, 0.25],2=>[0.2, 0.6, 0.2])
V = check_vc_increasing(fm)

gamma_ii = 1:5:10 #1:1:10
gamma_ij = 0.05:0.45:0.95 #0.05:0.1:0.95
a_j = linspace(12,300,2)

println(gamma_ii)
println(gamma_ij)
println(a_j)

elastic = false

SimulateOptimization(types, fm, a_j, gamma_ii, gamma_ij, elastic)

exit()







# # alphas
# a_1, a_2 = 0.25, 0.5
# # matrix for demand
# r_1, r_2 = 1,1.1
# gamma = 0.5
#
# boos = check_conditions_LM(types, fm, [a_1,a_2], [r_1,r_2], gamma)
#
#
# if r_1 + r_2 < 2*gamma
#     println("***ERROR: does not satisfy consistency check")
#     exit()
# end
#
# # build input matrices
# if length(fm) == 2
#     nalpha = [a_1; a_2]
#     nGamma = [r_1 -gamma;-gamma r_2]
# end
# if length(fm) == 1
#     nalpha = [a_1]
#     nGamma = [r_1]
# end
#
#


# version = "decentralized" # or decentralized
# elastic = true
#
# nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t = GenerateInputs(types, fm)
# D, wD, q_x, wq_x = InputsObjectiveLM(nalpha, nGamma, fm, nsupp, ntypes, nvars, sts)
# obj, transfers, x_vals, t_vals = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, nalpha, nGamma, fm, version, elastic)
#
# # Barrier solved model in 15 iterations and 0.02 seconds
# # Optimal objective 1.20162745e-01
# #
# # ===============================
# # Objective value: 0.12016274467115176
# # Allocations: [0.2, 0.44, 0.266667, 0.293333, 0.375, 0.0550002, 0.0500001, 0.515, 0.116667, 0.368333, 0.225, 0.13, 5.0723e-9, 0.54, 4.70934e-9, 0.426667, 0.0750012, 0.204999]
# # Transfers: [0.0679167, 0.124222, 0.0226389, 0.105778, 0.0679167, 0.0433333, 0.0429167, 0.0621111, 0.0143056, 0.0528889, 0.0429167, 0.0216667, 0.0062501, 0.124222, 0.00208337, 0.105778, 0.0062501, 0.0433333]
# # ===============================

# exit()
