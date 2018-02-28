include("utilities.jl")


function demand_LM(p, a, r, gamma, idx)
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

function asssortment_HM(p, loc, delta)
    out = []
    for i in 1:length(p)
        # if p[i] < min(p[j] + delta*abs(loc[j]-loc[i]) for j=1:length(p) if j!=i)
        aux = min([p[j] + delta*abs(loc[j]-loc[i]) for j=1:length(p) if j!=i]...)
        println(i, " ", aux)
        if p[i] < aux
            push!(out,i)
        end
    end
    return out

end

function demand_HM(p, loc, delta, idx, Q)
    # affine demand function for 2 suppliers
    if idx in Q
        if length(Q) == 1
            return 1
        end
        l_j, p_j = loc[idx], p[idx]
        if idx == min(Q...)
            idx_next = min([j for j in Q if j > idx]...)
            l_k, p_k = loc[idx_next], p[idx_next]
            aux_next = (p_k-p_j+delta*abs(l_k-l_j))/(2*delta)
            return l_j + aux_next
        elseif idx == max(Q...)
            idx_prev = max([j for j in Q if j < idx]...)
            l_i, p_i = loc[idx_prev], p[idx_prev]
            aux_prev = (p_i-p_j+delta*abs(l_j-l_i))/(2*delta)
            return (1-l_j) + aux_prev
        else
            idx_next = min([j for j in Q if j > idx]...)
            idx_prev = max([j for j in Q if j < idx]...)
            l_i, p_i = loc[idx_prev], p[idx_prev]
            l_k, p_k = loc[idx_next], p[idx_next]
            aux_next = (p_k-p_j+delta*abs(l_k-l_j))/(2*delta)
            aux_prev = (p_i-p_j+delta*abs(l_j-l_i))/(2*delta)
            return aux_prev + aux_next
        end
    else
        return 0
    end
    # p_i = 0
    # p_j = 0
    # if idx == 1
    #     p_i = p[1]
    #     p_j = p[2]
    # else
    #     p_i = p[2]
    #     p_j = p[1]
    # end
    # numerator = p_j-p_i+delta
    # denominator = 2*delta
    # return numerator/denominator
end




### INPUTS ###
# types = Dict(1=>0.1, 2=>0.2)
# fm = Dict(1=>[0.6,0.4],2=>[0.75,0.25])
types = Dict(1=>10, 2=>12)
fm = Dict(1=>[0.5,0.5],2=>[0.5,0.5])
# types = Dict(1=>0.1, 2=>0.2, 3=>0.25)
# fm = Dict(1=>[0.25,0.5, 0.25],2=>[0.2, 0.6, 0.2])

V = check_vc_increasing(fm)

### Linear Model Params ###
# alphas
a_1, a_2 = 0.25, 0.5
a = [a_1, a_2]

# matrix for demand
r_1, r_2 = 1,1.1
r = [r_1, r_2]
gamma = 0.5

### Hotelling Model Params ###
L = [0,1]
delta = 1


### MAIN ###
theta = [types[i] for i in keys(types)]
sort!(theta)
ntypes = length(types)
nsupp = length(fm)
Theta = combwithrep(nsupp, ntypes)


check_conditions_HM(Theta, fm, L, delta)
exit()


# compute cumulative distribution
Fm = ComputeCumulativeDistribution(fm)

total_cost = 0
for i in 1:length(Theta)
    prob = prod([fm[j][Theta[i][j]] for j=1:nsupp])
    p_s = []
    for j=1:nsupp
        idx = Theta[i][j]
        v_j = virtual_cost(theta, fm, Fm, j, idx)
        push!(p_s, v_j)
    end
    Q = asssortment_HM(p_s, L, delta)
    d_s = []
    for j=1:nsupp
        daux = demand_HM(p_s, L, delta, j, Q)
        push!(d_s, daux)
    end
    total_cost = total_cost + prob*d_s'*p_s
    #
    #
    # idx_s1, idx_s2 = Theta[i][1], Theta[i][2]
    # v_1 = virtual_cost(theta, fm, Fm, 1, idx_s1)
    # v_2 = virtual_cost(theta, fm, Fm, 2, idx_s2)
    # p_s = [v_1, v_2]
    # # d_1 = demand_LM(p_s, a, r, gamma, 1)
    # # d_2 = demand_LM(p_s, a, r, gamma, 2)
    # Q = asssortment_HM(p_s, L, delta)
    # println(Q)
    # d_1 = demand_HM(p_s, L, delta, 1, Q)
    # d_2 = demand_HM(p_s, L, delta, 2, Q)
    #
    # d_s = [d_1, d_2]
    println("========= Scenario ", i, " ==========")
    println("Types: ", Theta[i])
    println("Virtual costs: ", p_s)
    println("Demands: ", d_s)
end

println("TOTAL COST: ", total_cost)
