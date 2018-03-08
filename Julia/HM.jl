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

function CheckFeasibility(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, loc, delta, version, elastic=false, x0=nothing, t0=nothing)
    if version == "centralized"
        m = Model(solver=GurobiSolver(MIPGap = 1e-12))
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
    @variable(m, k[1:nvars])

    # constraints centralized version: IR, IC, feas
    @constraint(m, bG_x*x + bG_t*t .>= bh) # IR + IC
    if elastic == false
        @constraint(m, bA*x .== bb) # feas
    end
    @constraint(m, def_k[i=1:nsupp, j=1:sts],
                -0.5*delta*((loc[i] - sum(x[k] for k=((j-1)*nsupp+1):((j-1)*nsupp+i-1)))^2
                            + (sum(x[k] for k=((j-1)*nsupp+1):((j-1)*nsupp+i))-loc[i])^2 ) - k[(j-1)*nsupp+i]>= 0)

    if version == "decentralized"
        @variable(m, p[1:nvars]>=0)
        @variable(m, u[1:nvars]>=0)
        @variable(m, v[1:nvars])
        # @constraint(m, kkt_opt[i=1:nsupp, j=1:sts],
        #         p[(j-1)*nsupp+i] - u[(j-1)*nsupp+i] + v[(j-1)*nsupp+i]
        #         + delta*( sum( ( sum( x[s] for s=((j-1)*nsupp+1):r)-loc[r-(j-1)*nsupp] ) for r=((j-1)*nsupp+i):(j*nsupp) ) )
        #         - delta*( sum( ( loc[r+1-(j-1)*nsupp] - sum( x[s] for s=((j-1)*nsupp+1):r) ) for r=((j-1)*nsupp+i):(j*nsupp-1) ) ) == 0)
        @constraint(m, kkt_opt[i in 1:nvars], p[i] - u[i] + v[i] + delta*x[i] == 0)
        @NLconstraint(m, var_def[i in 1:nvars], t[i]-x[i]*p[i] == 0)
        @NLconstraint(m, kkt_comp[i in 1:nvars], x[i]*u[i] == 0)
        @constraint(m, kkt_cons[i in 1:sts], v[nsupp*(i-1)+1]-v[nsupp*(i-1)+2] == 0)
    end

    @objective(m, Max, 1)

    print(m)
    # exit()
    status = solve(m)
    obj = -getobjectivevalue(m)
    x_vals = getvalue(x)
    t_vals = getvalue(t)
    transfers = wq_t'*t_vals
    # println("Objective value: ", getobjectivevalue(m))
    #
    println("===============================")
    println("Objective value: ", getobjectivevalue(m))
    println("Allocations: ", getvalue(x))
    println("Transfers: ", getvalue(t))
    println("===============================")
    return obj, transfers, x_vals, t_vals
end

function SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, loc, delta, version, elastic=false, x0=nothing, t0=nothing)
    if version == "centralized"
        m = Model(solver=GurobiSolver(MIPGap = 1e-12))
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
        @variable(m, x[i=1:nvars]>=0, start=x0[i])
        @variable(m, t[i=1:nvars]>=0, start=t0[i])
    else
        @variable(m, x[1:nvars]>=0)
        @variable(m, t[1:nvars]>=0)
    end


    @variable(m, k[1:nvars])
    @variable(m, z)

    # constraints centralized version: IR, IC, feas
    @constraint(m, bG_x*x + bG_t*t .>= bh) # IR + IC
    if elastic == false
        @constraint(m, bA*x .== bb) # feas
    end
    @constraint(m, z <=  wq_t'*k - wq_t'*t)
    @constraint(m, def_k[i=1:nsupp, j=1:sts],
                -0.5*delta*((loc[i] - sum(x[k] for k=((j-1)*nsupp+1):((j-1)*nsupp+i-1)))^2
                            + (sum(x[k] for k=((j-1)*nsupp+1):((j-1)*nsupp+i))-loc[i])^2 ) - k[(j-1)*nsupp+i]>= 0)

    if version == "decentralized"
        @variable(m, p[1:nvars]>=0)
        @variable(m, u[1:nvars]>=0)
        if elastic == false
            @variable(m, v[1:nvars])
            @constraint(m, kkt_opt[i in 1:nvars], p[i] - u[i] + v[i] + delta*x[i] == 0)
            @constraint(m, kkt_cons[i in 1:sts], v[nsupp*(i-1)+1]-v[nsupp*(i-1)+2] == 0)
        else
            @constraint(m, kkt_opt[i in 1:nvars], p[i] - u[i] + delta*x[i] == 0)
        end
        @NLconstraint(m, var_def[i in 1:nvars], t[i]-x[i]*p[i] == 0)
        @NLconstraint(m, kkt_comp[i in 1:nvars], x[i]*u[i] == 0)

    end

    @objective(m, Max, z)

    status = solve(m)
    obj = -getobjectivevalue(m)
    x_vals = getvalue(x)
    t_vals = getvalue(t)
    transfers = wq_t'*t_vals
    # println("Objective value: ", getobjectivevalue(m))
    #
    println("===============================")
    println("Objective value: ", getobjectivevalue(m))
    println("Allocations: ", getvalue(x))
    println("Transfers: ", getvalue(t))
    println("===============================")
    return obj, transfers, x_vals, t_vals
end

function SimulateOptimization(types, fm, loc, deltas, elastic=false)
    println("Generating Inputs")
    nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t = GenerateInputs(types, fm)
    outfile = string("outputs/simulations_outcome_HM_", join(fm[1],'_'))
    if elastic
        outfile = string(outfile, "_elastic.txt")
    else
        outfile = string(outfile, "_inelastic.txt")
    end
    f = open(outfile, "w")
    for i=1:length(deltas)
        delta = deltas[i]
        println("Solving the problem for ", delta)
        obj_cent, tra_cent, x_cent, t_cent = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, loc, delta, "centralized", elastic)
        obj_dec, tra_dec, x_dec, t_dec = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, loc, delta, "decentralized", elastic, x_cent, t_cent)
        str_x_cent, str_t_cent = join(x_cent,';'), join(t_cent,';')
        str_x_dec, str_t_dec = join(x_dec,';'), join(t_dec,';')
        outstr = string(delta, ";", obj_cent, ";", str_x_cent, ";", str_t_cent,
                               ";", obj_dec, ";", str_x_dec, ";", str_t_dec)
        write(f, "$outstr \n")
    end
    close(f)
end


# IMPORTANT: types must be sorted in increasing order
### INPUTS ###
types = Dict(1=>10, 2=>12)
fm = Dict(1=>[0.6,0.4],2=>[0.6,0.4])
# types = Dict(1=>0.1, 2=>0.2)
# fm = Dict(1=>[0.6,0.4],2=>[0.75,0.25])
# types = Dict(1=>0.1, 2=>0.2, 3=>0.25)
# fm = Dict(1=>[0.25,0.5, 0.25],2=>[0.2, 0.6, 0.2])
V = check_vc_increasing(fm)

# alphas
loc = [0,1]
delta = 3
deltas = [i for i=0.5:0.5:6]

version = "decentralized" # or decentralized
elastic = false


elasticities = [false]
distributions = [[0.1, 0.9], [0.25, 0.75], [0.4, 0.6], [0.5,0.5], [0.6,0.4], [0.75, 0.25], [0.9, 0.1]]

# nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t = GenerateInputs(types, fm)
# obj_cent, tra_cent, x_cent, t_cent = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, loc, delta, "centralized", elastic)
# obj_dec, tra_dec, x_dec, t_dec = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, loc, delta, "decentralized", elastic, x_cent, t_cent)

for e=1:length(elasticities)
    elastic = elasticities[e]
    for d=1:length(distributions)
        distr = distributions[d]
        fm = Dict(1=>distr,2=>distr)
        SimulateOptimization(types, fm, loc, deltas, elastic)
    end
end


exit()

nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t = GenerateInputs(types, fm)
obj_cent, tra_cent, x_cent, t_cent = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, loc, delta, "centralized")




# Objective value: -12.125009520159871
# Allocations: [0.499992, 0.500008, 0.999999, 7.4766e-7, 5.79408e-7, 0.999999, 0.500003, 0.499997]
# Transfers: [15.9999, 16.0001, 0.0, 5.99998, 6.00004, 0.0, 0.0, 0.0]


# for i=1:length(u0)
#     println(p0[i] - u0[i] + v0[i] + delta*x_cent[i] )
# end

# obj_cent, tra_cent, x_dec, t_dec = CheckFeasibility(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, loc, delta, "decentralized", x_cent, t_cent)

# obj_cent, tra_cent, x_cent, t_cent = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, loc, delta, "decentralized", x_cent, t_cent)

# Objective value: -12.458331756709534
# Allocations: [0.833335, 0.166665, 0.833335, 0.166665, 0.0, 1.0, 0.0, 1.0]
# Transfers: [8.33334, 1.99998, 8.33335, 1.99998, 0.0, 12.0, 0.0, 12.0]

exit()
