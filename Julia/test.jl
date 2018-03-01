# using JuMP
# using Gurobi, KNITRO
include("utilities.jl")

using PyPlot
# x = linspace(0,2*pi,1000); y = sin.(3*x + 4*cos.(2*x))
# plot(x, y, color="red", linewidth=2.0, linestyle="--")
# show()


### -------------------------- ###
### INPUTS ###
### -------------------------- ###
# types = Dict(1=>0.1, 2=>0.2)
# fm = Dict(1=>[0.6,0.4],2=>[0.75,0.25])
types = Dict(1=>10, 2=>12)
fm = Dict(1=>[0.5,0.5],2=>[0.5,0.5])
# types = Dict(1=>0.1, 2=>0.2, 3=>0.25)
# fm = Dict(1=>[0.25,0.5, 0.25],2=>[0.2, 0.6, 0.2])

model_to_test = "HM" # "HM" for Hotelling mdel or "LM" for linear model

### Hotelling Model Params ###
L = [0,1]
deltas = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7]


theta = [types[i] for i in keys(types)]
sort!(theta)
ntypes = length(types)
nsupp = length(fm)
Theta = combwithrep(nsupp, ntypes)

# check for assumptions and conditions
V = check_vc_increasing(fm)
Fm = ComputeCumulativeDistribution(fm)





costs = []
for s=1:length(deltas)
    delta = deltas[s]
    boos = check_conditions_HM(Theta, fm, L, delta)

    total_cost = 0
    for i in 1:length(Theta)
        prob = prod([fm[j][Theta[i][j]] for j=1:nsupp])
        p_s = [V[j,Theta[i][j]] for j=1:nsupp]
        Q = asssortment_HM(p_s, L, delta)
        d_s = []
        for j=1:nsupp
            daux = demand_HM(p_s, L, delta, j, Q)
            push!(d_s, daux)
        end
        total_cost = total_cost + prob*d_s'*p_s
    end
    push!(costs, total_cost)
end

plot(deltas, costs, color="red", linewidth=2.0, linestyle="--")
show()
