using JuMP
using PyPlot
using Gurobi, KNITRO
include("utilities.jl")
include("LDM.jl")

# ----------------
# SIMULATION METHODS
# ----------------
# Two suppliers, asymmetric (i.e. parameters can be different)
function SimulateTwoSuppliersAsymmetric(types, fm, a_j, gamma_ii, gamma_ij, elastic=false, expostir=false)
    ## OVERVIEW ##
    # This method simulates a two supplier general affine model for different combinations of
    # values of qualities and demand matrices. It considers assymetric settings, i.e.
    # cases where suppliers qualities are different and when the demand matrix is not
    # symmetric.

    ## INPUTS ##
    # 1. types = dictionary with types, enumerated from 1 to n
    # 2. fm = dictionary that tells, for each supplier, what is the distribution of types
    # 3. a_j = list of possible qualities we want to simulate
    # 4. gamma_ii = list of possible values in the diagonal of demand matrix
    # 5. gamma_ij = list of possible values in the off-diagonal of demand matrix
    # 6. elastic = boolean; true => elastic demand, false=> inelastic demand
    # 7. expostir = boolean; true => expostir constraints, false=> only interim constraints

    # Check whether inputs are correct
    if length(fm) != N
        println("***ERROR: the number of suppliers does not match with the number of marginal distributions")
        println("Number of suppliers:", N)
        println("Numver of marginal distributions: ", length(fm))
        exit()
    end

    if length(types) != length(fm[1])
        println("***ERROR: the number of types does not match the length of the distribution list")
        println("Number of types:", length(types))
        println("Numver of elements in distribution: ", length(dist))
        exit()
    end
    V = check_vc_increasing(fm)

    outfile = string("outputs/simulations_outcome_LM_asymmetric_", join(fm[1],'_'))
    if expostir
        outfile = string(outfile, "_expostir")
    end
    if elastic
        outfile = string(outfile, "_elastic.txt")
    else
        outfile = string(outfile, "_inelastic.txt")
    end
    f = open(outfile, "w")
    for i1=1:length(a_j)
        # to start assume that qualities are the same for both suppliers and same price sensitivities
        for i2=i1:length(a_j)
            alpha = [a_j[i1]; a_j[i2]]
            for j1=1:length(gamma_ii)
                for j2=j1:length(gamma_ii)
                    for k=1:length(gamma_ij)
                        Gamma = [gamma_ii[j1] -gamma_ij[k]; -gamma_ij[k] gamma_ii[j2]]
                        if prod([check_diagonal_dominant(Gamma), check_consistency_demand_matrix(Gamma), isposdef(Gamma)]) == false
                            println("***WARNING: the matrix does not satisfy the assumptions")
                            continue
                        end
                        obj_cent, transfers_cent, x_cent, t_cent, obj_dec, transfers_dec, x_dec, t_dec = FormulateAndSolve(types, fm, alpha, Gamma, elastic, expostir)

                        str_x_cent, str_t_cent = join(x_cent,';'), join(t_cent,';')
                        str_x_dec, str_t_dec = join(x_dec,';'), join(t_dec,';')

                        outstr = string(a_j[i1], ";", a_j[i2], ";", gamma_ii[j1], ";", gamma_ii[j2], ";", gamma_ij[k],
                                               ";", obj_cent, ";", str_x_cent, ";", str_t_cent,
                                               ";", obj_dec, ";", str_x_dec, ";", str_t_dec)
                        write(f, "$outstr \n")
                    end
                end
            end
        end
    end
    close(f)
end

# N suppliers
function SimulateNSupplierSymmetric(types, fm, a_j, gamma_ii, gamma_ij, elastic=false, expostir=false, N=2)
    ## OVERVIEW ##
    # This method simulates a N supplier general affine model for different combinations of
    # values of qualities and demand matrices. It considers symmetric settings, i.e.
    # cases where suppliers qualities are the same and the demand matrix is
    # symmetric.

    ## INPUTS ##
    # 1. types = dictionary with types, enumerated from 1 to n
    # 2. fm = dictionary that tells, for each supplier, what is the distribution of types
    # 3. a_j = list of possible qualities we want to simulate
    # 4. gamma_ii = list of possible values in the diagonal of demand matrix
    # 5. gamma_ij = list of possible values in the off-diagonal of demand matrix
    # 6. elastic = boolean; true => elastic demand, false=> inelastic demand
    # 7. expostir = boolean; true => expostir constraints, false=> only interim constraints

    # Check whether inputs are correct
    if length(fm) != N
        println("***ERROR: the number of suppliers does not match with the number of marginal distributions")
        println("Number of suppliers:", N)
        println("Numver of marginal distributions: ", length(fm))
        exit()
    end

    if length(types) != length(fm[1])
        println("***ERROR: the number of types does not match the length of the distribution list")
        println("Number of types:", length(types))
        println("Numver of elements in distribution: ", length(dist))
        exit()
    end
    V = check_vc_increasing(fm)


    # start process
    outfile = string("outputs/simulations_outcome_LM_symmetric_N=", N,"_", join(fm[1],'_'))
    if expostir
        outfile = string(outfile, "_expostir")
    end
    if elastic
        outfile = string(outfile, "_elastic.txt")
    else
        outfile = string(outfile, "_inelastic.txt")
    end
    f = open(outfile, "w")
    for i1=1:length(a_j)
        alpha = a_j[i1]*ones(N)
        for j1=1:length(gamma_ii)
            for k=1:length(gamma_ij)
                Gamma =  gamma_ii[j1]*eye(N)-gamma_ij[k]*(ones(N,N)-eye(N))
                if prod([check_diagonal_dominant(Gamma), check_consistency_demand_matrix(Gamma), isposdef(Gamma)]) == false
                    println("***WARNING: the matrix does not satisfy the assumptions")
                    continue
                end
                obj_cent, transfers_cent, x_cent, t_cent, obj_dec, transfers_dec, x_dec, t_dec = FormulateAndSolve(types, fm, alpha, Gamma, elastic, expostir)

                str_x_cent, str_t_cent = join(x_cent,';'), join(t_cent,';')
                str_x_dec, str_t_dec = join(x_dec,';'), join(t_dec,';')

                outstr = string(a_j[i1], ";", a_j[i1], ";", gamma_ii[j1], ";", gamma_ii[j1], ";", gamma_ij[k],
                                       ";", obj_cent, ";", str_x_cent, ";", str_t_cent,
                                       ";", obj_dec, ";", str_x_dec, ";", str_t_dec)
                write(f, "$outstr \n")
            end
        end
    end
    close(f)
end

# =================
# EXECUTION
# =================
# TO SET PARAMETERS FOR SIMULATION, CHANGE HERE!

# SUPPLIERS AND TYPES
N = 2 # number of suppliers
types = Dict(1=>10, 2=>10.5)
dist = [0.2,0.8]
fm = Dict(i=>dist for i=1:N)


# QUALITIES AND DEMAND PARAMETERS
# if you want to run a single instance, enter a list with just one element
gamma_ii = [1] #1:1:10 #4
gamma_ij = [0.05]#0.05:0.1:0.95 #0.05:0.1:0.95 #0.05:0.1:0.95 #4
a_j = [11]# linspace(25,100,4)

# SETTING FOR SIMULATIONS
elasticities = [false, true] #, false
expostirs = [false, true]
distributions = [[0.1, 0.9], [0.2, 0.8], [0.25, 0.75], [0.4, 0.6], [0.5, 0.5], [0.6,0.4], [0.75, 0.25], [0.8, 0.2], [0.9,0.1]]

for e=1:length(elasticities)
    elastic = elasticities[e]
    for d=1:length(distributions)
        distr = distributions[d]
        fm = Dict(i=>distr for i=1:N) # two suppliers
        for ex=1:length(expostirs)
            expostir = expostirs[ex]
            SimulateNSupplierSymmetric(types, fm, a_j, gamma_ii, gamma_ij, elastic, expostir, N)
        end
    end
end
