using JuMP
using PyPlot
using Gurobi, KNITRO
include("utilities.jl")
include("HM.jl")

# ----------------
# SIMULATION METHODS
# ----------------
# Two suppliers, asymmetric (i.e. parameters can be different)
function SimulateTwoSuppliers(types, fm, loc, deltas, elastic=false, expostir=false, outdir="outputs")
    if !isdir(outdir)
        mkdir(outdir)
    end

    println("Generating Inputs")
    nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, fth, Theta = GenerateInputs(types, fm)
    outfile = string(outdir,"/simulations_outcome_HM_", join(fm[1],'_'))
    if expostir
        outfile = string(outfile, "_expostir")
    end
    if elastic
        outfile = string(outfile, "_elastic.txt")
    else
        outfile = string(outfile, "_inelastic.txt")
    end
    println(outfile)
    f = open(outfile, "w")
    out = []
    for i=1:length(deltas)
        delta = deltas[i]
        println("Solving the problem for ", delta)
        obj_cent, tra_cent, x_cent, t_cent = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, types, Theta, loc, delta, "centralized", elastic, expostir)
        obj_dec, tra_dec, x_dec, t_dec = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, types, Theta, loc, delta, "decentralized", elastic, expostir, x_cent, t_cent)
        str_x_cent, str_t_cent = join(x_cent,';'), join(t_cent,';')
        str_x_dec, str_t_dec = join(x_dec,';'), join(t_dec,';')
        outstr = string(delta, ";", obj_cent, ";", str_x_cent, ";", str_t_cent,
                               ";", obj_dec, ";", str_x_dec, ";", str_t_dec)
        write(f, "$outstr \n")
        push!(out, "$outstr \n")
    end
    close(f)
    for i=1:length(out)
        println(out[i])
    end
end

# =================
# EXECUTION
# =================
# TO SET PARAMETERS FOR SIMULATION, CHANGE HERE!
# SUPPLIERS AND TYPES
# types = Dict(1=>8, 2=>10, 3=>12)
# fm = Dict(1=>[0.25,0.5,0.25],2=>[0.25,0.5,0.25])

types = Dict(1=>10, 2=>12)
fm = Dict(1=>[0.25,0.75],2=>[0.25,0.75])

# types = Dict(1=>0.1, 2=>0.2)
# fm = Dict(1=>[0.6,0.4],2=>[0.75,0.25])
# types = Dict(1=>0.1, 2=>0.2, 3=>0.25)
# fm = Dict(1=>[0.25,0.5, 0.25],2=>[0.2, 0.6, 0.2])
V = check_vc_increasing(fm)

# PARAMETERS OF THE PROBLEM: LOCATIONS AND TRANSPORTATION COST
loc = [0,1]
delta = 0.1

nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, fth, Theta = GenerateInputs(types, fm)
obj_cent, tra_cent, x_cent, t_cent = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, types, Theta, loc, delta, "centralized", false, false)


exit()


elasticities = [false] #, false
expostirs = [false, true]
distributions = [[(1/3), (1/3), (1/3)],
                 [0.25,0.25, 0.5], [0.25, 0.5, 0.25], [0.5, 0.25, 0.25],
                 [0.2, 0.2, 0.6], [0.2, 0.6, 0.2], [0.6, 0.2, 0.2],
                 [0.1, 0.1, 0.8], [0.1, 0.8, 0.1], [0.8, 0.1, 0.1]]

for e=1:length(elasticities)
    elastic = elasticities[e]
    for d=1:length(distributions)
        distr = distributions[d]
        fm = Dict(1=>distr,2=>distr) # two suppliers
        for ex=1:length(expostirs)
            expostir = expostirs[ex]
            SimulateTwoSuppliers(types, fm, loc, deltas, elastic, expostir)
        end
    end
end
