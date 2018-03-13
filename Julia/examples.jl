include("utilities.jl")

### -------------------------- ###
### INPUTS ###
### -------------------------- ###
# types = Dict(1=>0.1, 2=>0.2)
# fm = Dict(1=>[0.6,0.4],2=>[0.75,0.25])
# types = Dict(1=>10, 2=>12)
# fm = Dict(1=>[0.5,0.5],2=>[0.5,0.5])
# types = Dict(1=>0.1, 2=>0.2, 3=>0.25)
# fm = Dict(1=>[0.25,0.5, 0.25],2=>[0.2, 0.6, 0.2])

types = Dict(1=>10, 2=>12)
# fm = Dict(1=>[0.1,0.9],2=>[0.1,0.9])
fm = Dict(1=>[0.5,0.5],2=>[0.5,0.5])

model_to_test = "HM" # "HM" for Hotelling mdel or "LM" for linear model

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
delta = 2

### -------------------------- ###
### MAIN ###
### -------------------------- ###
theta = [types[i] for i in keys(types)]
sort!(theta)
ntypes = length(types)
nsupp = length(fm)
Theta = combwithrep(nsupp, ntypes)

# check for assumptions and conditions
V = check_vc_increasing(fm)




# compute cumulative distribution
Fm = ComputeCumulativeDistribution(fm)

total_cost = 0
total_objective = 0
if model_to_test == "HM"
    boos = check_conditions_HM(Theta, fm, L, delta)
    for i in 1:length(Theta)
        println("========= Scenario ", i, " ==========")
        prob = prod([fm[j][Theta[i][j]] for j=1:nsupp])
        p_s = [V[j,Theta[i][j]] for j=1:nsupp]

        Q = asssortment_HM(p_s, L, delta)
        d_s = []
        for j=1:nsupp
            daux = demand_HM(p_s, L, delta, j, Q)
            push!(d_s, daux)
        end
        total_cost = total_cost + prob*d_s'*p_s
        total_objective = total_objective + prob*d_s'*p_s + 0.5*delta*prob*sum([d_s[j]^2 for j=1:nsupp])
        println("Types: ", Theta[i])
        println("Probability: ", prob)
        println("Virtual costs: ", p_s)
        println("Demands: ", d_s)
    end
elseif model_to_test == "LM"
    boos = check_conditions_LM(Theta, fm, a, r, gamma)
    for i in 1:length(Theta)
        println("========= Scenario ", i, " ==========")
        prob = prod([fm[j][Theta[i][j]] for j=1:nsupp])
        p_s = [V[j,Theta[i][j]] for j=1:nsupp]
        d_s = []
        for j=1:nsupp
            daux = demand_LM(p_s, a, r, gamma, j)
            push!(d_s, daux)
        end
        total_cost = total_cost + prob*d_s'*p_s
        println("Types: ", Theta[i])
        println("Virtual costs: ", p_s)
        println("Demands: ", d_s)
    end
else
    println("***ERROR: unknown model to test. Please consider either LM or HM.")
end

println("")
println("=================================")
println("TOTAL OBJECTIVE: ", total_objective)
println("TOTAL TRANSFERS: ", total_cost)
println("=================================")
