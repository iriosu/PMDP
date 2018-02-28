using Iterators

types = Dict(1=>0.1, 2=>0.2)
# println(types)
# for i in keys(types)
#     println(i, types[i])
# end
ntypes = length(types)
fm = Dict(1=>[0.6,0.4],2=>[0.75,0.25])
nsupp = length(fm)
# println(nsupp, ntypes)
Theta = []
for p in product(1:ntypes,1:ntypes)
    push!(Theta,p)
end
# println(Theta)
# println(Theta[end][1])
# joint distribution based on marginals
f = Dict()
for perm in Theta
    f[perm] = prod([fm[i][perm[i]] for i=1:nsupp])
end

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
if nsupp == 2
    nalpha = [a_1; a_2]
    nGamma = [r_1 -gamma;-gamma r_2]
end
if nsupp == 1
    nalpha = [a_1]
    nGamma = [r_1]
end
nD = inv(nGamma)
# println(size(nD))
# println(size(nalpha))

nc = nD*nalpha
# println(nc)
# println(size(nc))


# ------------------------------------------------------------------------------
# CENTRALIZED WITH IC AND IR
# ------------------------------------------------------------------------------
nvars = length(Theta)*nsupp # here we include transfers and allocations
sts = length(Theta) # size of type space
x0 = ones(nvars)
t0 = ones(nvars)
# println(x0)


# this dict tell us the probability of all other subjects being of their type $f_{-i}(\theta_{-i})$
f_woi = Dict(i=>Dict(j=>0.0 for j in Theta) for i in 1:nsupp)
# println(f_woi)
# println(Theta)


for i in keys(f_woi)
    println(i)
    supp_woi = [j for j in 1:nsupp if j!=i]

    # println(supp_woi)
    for perm in Theta
        f_woi[i][perm] = prod([fm[j][perm[i]] for j in supp_woi])
    end
end
# println(f_woi)
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
        println(i," ", j," ", idxs)
        for k in idxs
            bIR_x[row, nsupp*(k-1)+i] = -types[j]*f_woi[i][Theta[k]]
            bIR_t[row, nsupp*(k-1)+i] = f_woi[i][Theta[k]]
        end
        row=row+1
    end
end
println(bIR_x)
println(bIR_t)


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
# println(bG)
println(size(bG_x))
println(size(bG_t))

bD = zeros((nvars, nvars))
bA = zeros((length(Theta), nvars))
bh = zeros(size(bG_x)[1])
bb = ones(length(Theta))

bq_x = zeros(nvars)
bq_t = zeros(nvars)


for i in 1:length(Theta)
    bD[nsupp*(i-1)+1:nsupp*i,nsupp*(i-1)+1:nsupp*i] = f[Theta[i]]*nD
    bA[i,nsupp*(i-1)+1:nsupp*i] = ones(nsupp)
    bq_x[nsupp*(i-1)+1:nsupp*i] = f[Theta[i]]*nc
    bq_t[nsupp*(i-1)+1:nsupp*i] = -f[Theta[i]]*ones(nsupp)
end

println(bA)

aux = 0.5*x0'*bD*x0
aux2 = bq_x'*x0
println(aux2)
I = eye(nvars)
e = ones(nvars)




# OPTIMIZATION
using JuMP
using Gurobi, KNITRO


# m = Model(solver=GurobiSolver(Presolve=0))
m = Model(solver=KnitroSolver(mip_method = KTR_MIP_METHOD_BB,
                              ms_enable = 1, ms_maxsolves = 1000,
                              algorithm = KTR_ALG_ACT_CG,
                              outmode = KTR_OUTMODE_SCREEN,
                              KTR_PARAM_OUTLEV = KTR_OUTLEV_ALL,
                              KTR_PARAM_MIP_OUTINTERVAL = 1,
                              KTR_PARAM_MIP_MAXNODES = 10000,
                              KTR_PARAM_HESSIAN_NO_F = KTR_HESSIAN_NO_F_ALLOW))
@variable(m, x[1:nvars]>=0)
@variable(m, p[1:nvars]>=0)
@variable(m, t[1:nvars]>=0)
@variable(m, u[1:nvars]>=0)
@variable(m, v[1:nvars])
@variable(m, z)
@constraint(m, bG_x*x + bG_t*t .<= bh) # IR + IC
@constraint(m, bA*x .== bb) # feas
@constraint(m, kkt_opt, p - bq_x + bD*x - u + v .== 0)
@NLconstraint(m, var_def[i in 1:nvars], t[i]-x[i]*p[i] == 0)
@NLconstraint(m, kkt_comp[i in 1:nvars], x[i]*u[i] == 0)
@constraint(m, kkt_cons[i in 1:sts], v[nsupp*(i-1)+1]-v[nsupp*(i-1)+2] == 0)
@constraint(m, z <=  bq_x'*x + bq_t'*t - 0.5*x'*bD*x )
@objective(m, Max, z)

print(m)

status = solve(m)

println("Objective value: ", getobjectivevalue(m))

exit()
