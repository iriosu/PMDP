include("utilities.jl")

# demand function in 2x2 example

function demand(p, a, r, gamma, idx)
    r_i, a_i, p_i = 0, 0, 0
    r_j, a_j, p_j = 0, 0, 0
    if idx == 1
        r_i, a_i, p_i = r[1], a[1], p[1]
        r_j, a_j, p_j = r[2], a[2], p[2]
    else
        r_i, a_i, p_i = r[2], a[2], p[2]
        r_j, a_j, p_j = r[1], a[1], p[1]
    end
    numerator = (r_j-gamma)*a_i - (r_i-gamma)*a_j + r_i - gamma - (r_i*r_j-gamma^2)*(p_i-p_j)
    denominator = r_i + r_j - 2*gamma
    return max(0, min( numerator/denominator, 1 )  )

end


function virtual_cost(theta, fm, Fm, idx_s, idx_t)
    theta_i = theta[idx_t]
    f_i, F_i = fm[idx_s], Fm[idx_s]
    if idx_t == 1
        rho_i = 1
    else
        rho_i = idx_t - 1
    end
    # println("    ", theta_i, " ", F_i[rho_i], " ", f_i[idx_t], " ", theta[rho_i], " ", rho_i)
    return theta_i + (F_i[rho_i]/f_i[idx_t])*(theta_i - theta[rho_i])
end



### INPUTS ###
# types = Dict(1=>0.1, 2=>0.2)
# fm = Dict(1=>[0.6,0.4],2=>[0.75,0.25])

types = Dict(1=>0.1, 2=>0.2, 3=>0.25)
fm = Dict(1=>[0.5,0.25, 0.25],2=>[0.6,0.2, 0.2])

# alphas
a_1, a_2 = 0.25, 0.5
a = [a_1, a_2]

# matrix for demand
r_1, r_2 = 1,1.1
r = [r_1, r_2]

gamma = 0.5



### MAIN ###
theta = [types[i] for i in keys(types)]
sort!(theta)
ntypes = length(types)
nsupp = length(fm)
Theta = combwithrep(nsupp, ntypes)


# compute cumulative marginal distributions
Fm = Dict(i=>[] for i in keys(fm))
for key in keys(fm)
    for i in 1:length(fm[key])
        push!(Fm[key],sum([fm[key][j] for j in 1:i]))
    end
end

println(Fm)
println(theta)

for i in 1:nsupp
    for j in 1:ntypes
        println("Supplier ", i, " Type ", theta[j])
        println(virtual_cost(theta, fm, Fm, i, j))
    end
end

for i in 1:length(Theta)
    idx_s1, idx_s2 = Theta[i][1], Theta[i][2]
    v_1 = virtual_cost(theta, fm, Fm, 1, idx_s1)
    v_2 = virtual_cost(theta, fm, Fm, 2, idx_s2)
    p_s = [v_1, v_2]
    d_1 = demand(p_s, a, r, gamma, 1)
    d_2 = demand(p_s, a, r, gamma, 2)
    d_s = [d_1, d_2]
    println("========= Scenario ", i, " ==========")
    println("Types: ", Theta[i])
    println("Virtual costs: ", p_s)
    println("Demands: ", d_s)
end
