include("utilities.jl")


function demand(p, a, r, gamma, idx)
    # affine demand function for 2 suppliers
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



### INPUTS ###
# types = Dict(1=>0.1, 2=>0.2)
# fm = Dict(1=>[0.6,0.4],2=>[0.75,0.25])

types = Dict(1=>0.1, 2=>0.2, 3=>0.25)
fm = Dict(1=>[0.25,0.5, 0.25],2=>[0.2, 0.6, 0.2])
V = check_vc_increasing(fm)

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

# compute cumulative distribution
Fm = ComputeCumulativeDistribution(fm)

# we check if the assumptions of the model are satisfied, i.e. that v(theta) is increasing



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
