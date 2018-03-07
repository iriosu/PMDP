# ==============================
# Formulates the problem and solves it
# ==============================

include("utilities.jl")

function check_diagonal_dominant(A)
    boos = 2*diag(abs.(A)).>sum(abs.(A),2)
    if prod(boos) == false
        println("***WARNING: the matrix provided is not strictly diagonal dominant")
    end
    return prod(boos)
end

function check_consistency_demand_matrix(A)
    if sum(diag(A)) < 2*A[1,2]
        println("***ERROR: does not satisfy consistency check")
        return false
    else
        return true
    end
end

gamma_ii = 1:5:10 #1:1:10
gamma_ij = 0.05:0.45:0.95 #0.05:0.1:0.95
a_j = linspace(12,300,2)

println(gamma_ii)
println(gamma_ij)
println(a_j)

elastic = false

#nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t = GenerateInputs(types, fm)

if elastic
    f = open("simulations_outcome_elastic.txt", "w")
else
    f = open("simulations_outcome_inelastic.txt", "w")
end
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
                    outstr = string(a_j[i1], ";", a_j[i2], ";", gamma_ii[j1], ";", gamma_ii[j2], ";", gamma_ij[k])
                    write(f, "$outstr \n")
                    #obj_cent, transfers_cent, x_cent, t_cent = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, alpha, Gamma, "centralized", elastic)
                    #obj_dec, transfers_dec, x_dec, t_dec = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, alpha, Gamma, "decentralized", elastic, x_cent, t_cent)
                end
            end
        end
    end
end
close(f)
