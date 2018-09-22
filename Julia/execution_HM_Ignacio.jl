using DataFrames, Ipopt, JuMP, NLopt
using DataFrames, CSV
include("HM_Nlopt.jl")
directory = "/Users/ziminzeng/Dropbox/2018_Spring/RA/code_v1_HM/Bin/"

function three_types_original()
   loc = [0,1]
   # Parameters
   # types = Dict(1=>10, 2=>12)
   elastic = false
   types = Dict(1=>10, 2=>12, 3=>15)
   distributions = [[0.1, 0.2, 0.7],
                   [0.25, 0.35, 0.4],
                   [0.5, 0.4, 0.1],
                   [0.75, 0.2, 0.05],
                   [0.9, 0.035, 0.065]]
   deltas = collect(0.1:1:10)
   output = DataFrame(theta1 = Float64[],
                                    theta2 = Float64[],
                                    theta3 = Float64[],
                                    f1 = Float64[],
                                    f2 = Float64[],
                                    f3 = Float64[],
                                    delta = Float64[],
                                    cent = Float64[],
                                    decent = Float64[],
                                    decent_expostIR = Float64[])

   for delta in deltas
      for d=1:length(distributions)
         distr = distributions[d]
         fm = Dict(1=>distr,2=>distr) # two suppliers
         nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, fth, Theta = GenerateInputs(types, fm)
         cent, tra_cent, x_cent, t_cent = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, types, Theta, loc, delta, "centralized", false, false)
         decent, tra_dec, x_dec, t_dec = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, types, Theta, loc, delta, "decentralized", false, false, x_cent, t_cent)
         expostIR, tra_dec, x_dec, t_dec = SolveOptimization(nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, types, Theta, loc, delta, "decentralized", false , true, x_cent, t_cent)
         push!(output, [10, 12, 15, distr[1], distr[2], distr[3], delta, cent, decent, expostIR])
      end
   end
   filename = "ignacio.csv"
   CSV.write(directory*filename, output)
end

three_types_original()
