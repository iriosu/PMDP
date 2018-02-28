using JuMP
using Gurobi, KNITRO
include("utilities.jl")

# ==============================
# Formulates the problem and solves it
# ==============================
function FormulateAndSolve(types, fm, nalpha, nGamma, version)

    ntypes = length(types)
    nsupp = length(fm)
    Theta = combwithrep(nsupp, ntypes)

    # joint distribution based on marginals
    f = Dict()
    for perm in Theta
        f[perm] = prod([fm[i][perm[i]] for i=1:nsupp])
    end

    nD = inv(nGamma)
    # println(size(nD))
    # println(size(nalpha))

    nc = nD*nalpha
    # println(nc)
    # println(size(nc))

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
