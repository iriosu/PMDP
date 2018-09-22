using JuMP
# using PyPlotPkg.add("Clp")
using Gurobi, KNITRO, Ipopt, NLopt

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
    # CENTRALIZEPkg.add("Clp")D WITH IC AND IR
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

function SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, types, Theta, loc, delta, version, elastic=false, expostir=false, x0=nothing, t0=nothing)
    if version == "centralized"
        m = Model(solver=GurobiSolver(MIPGap = 1e-12))
    elseif version == "decentralized"
        m = Model(solver=NLoptSolver(algorithm=:LD_MMA))
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

    # normal constraints
    @constraint(m, z <=  wq_t'*k - wq_t'*t)
    @constraint(m, def_k[i=1:nsupp, j=1:sts],
                -0.5*delta*((loc[i] - sum(x[k] for k=((j-1)*nsupp+1):((j-1)*nsupp+i-1)))^2
                            + (sum(x[k] for k=((j-1)*nsupp+1):((j-1)*nsupp+i))-loc[i])^2 ) - k[(j-1)*nsupp+i]>= 0)
    @constraint(m, bG_x*x + bG_t*t .>= bh) # IR + IC interim

    # special constraints
    if expostir == true
        epG_x, epG_t, eph = GenerateExPostInputs(types, Theta, nsupp, nvars)
        @constraint(m, epG_x*x + epG_t*t .>= eph) # IR + IC
    end
    if elastic == false
        @constraint(m, bA*x .>= bb) # feas
        @constraint(m, bA*x .<= bb) # feas
    end


    if version == "decentralized"
        @variable(m, p[1:nvars]>=0)
        @variable(m, u[1:nvars]>=0)
        if elastic == false
            @variable(m, v[1:nvars])
            @constraint(m, kkt_opt[i in 1:nvars], p[i] - u[i] + v[i] + delta*x[i] >= 0)
            @constraint(m, kkt_opt[i in 1:nvars], p[i] - u[i] + v[i] + delta*x[i] <= 0)
            @constraint(m, kkt_cons[i in 1:sts], v[nsupp*(i-1)+1]-v[nsupp*(i-1)+2] >= 0)
            @constraint(m, kkt_cons[i in 1:sts], v[nsupp*(i-1)+1]-v[nsupp*(i-1)+2] <= 0)
        else
            @constraint(m, kkt_opt[i in 1:nvars], p[i] - u[i] + delta*x[i] >= 0)
            @constraint(m, kkt_opt[i in 1:nvars], p[i] - u[i] + delta*x[i] <= 0)
        end
        @NLconstraint(m, var_def[i in 1:nvars], t[i]-x[i]*p[i] >= 0)
        @NLconstraint(m, var_def[i in 1:nvars], t[i]-x[i]*p[i] <= 0)
        @NLconstraint(m, kkt_comp[i in 1:nvars], x[i]*u[i] >= 0)
        @NLconstraint(m, kkt_comp[i in 1:nvars], x[i]*u[i] <= 0)

    end

    @objective(m, Max, z)

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




# =================
# EXECUTION
# =================
# TO SET PARAMETERS FOR SIMULATION, CHANGE HERE!
# SUPPLIERS AND TYPES
# types = Dict(1=>8, 2=>10, 3=>12)
# fm = Dict(1=>[0.25,0.5,0.25],2=>[0.25,0.5,0.25])
types = Dict(1=>8, 2=>10, 3=>12)
fm = Dict(1=>[1/3,1/3,1/3],2=>[1/3,1/3,1/3])
# types = Dict(1=>0.1, 2=>0.2)
# fm = Dict(1=>[0.6,0.4],2=>[0.6,0.4])
# types = Dict(1=>0.1, 2=>0.2, 3=>0.25)
# fm = Dict(1=>[0.25,0.5, 0.25],2=>[0.2, 0.6, 0.2])
V = check_vc_increasing(fm)

# PARAMETERS OF THE PROBLEM: LOCATIONS AND TRANSPORTATION COST
loc = [0,1]
delta = 5
elastic = false
expostir = false


nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, fth, Theta = GenerateInputs(types, fm)
obj_cent, tra_cent, x_cent, t_cent = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, types, Theta, loc, delta, "centralized", elastic, expostir)
obj_dec, tra_dec, x_dec, t_dec = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, types, Theta, loc, delta, "decentralized", elastic, expostir, x_cent, t_cent)
