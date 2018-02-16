using JuMP


println("Starting")
# m = Model(solver = ClpSolver())
# @variable(m, 0 <= x <= 2 )
# @variable(m, 0 <= y <= 30 )
#
# @objective(m, Max, 5x + 3*y )
# @constraint(m, 1x + 5y <= 3.0 )
#
# print(m)
#
# status = solve(m)
#
# println("Objective value: ", getobjectivevalue(m))
# println("x = ", getvalue(x))
# println("y = ", getvalue(y))


using KNITRO, JuMP, Ipopt

## Solve test problem 1 (Synthesis of processing system) in
 #  M. Duran & I.E. Grossmann, "An outer approximation algorithm for
 #  a class of mixed integer nonlinear programs", Mathematical
 #  Programming 36, pp. 307-339, 1986.  The problem also appears as
 #  problem synthes1 in the MacMINLP test set.

m = Model(solver=KnitroSolver(mip_method = KTR_MIP_METHOD_BB,
                              algorithm = KTR_ALG_ACT_CG,
                              outmode = KTR_OUTMODE_SCREEN,
                              KTR_PARAM_OUTLEV = KTR_OUTLEV_ALL,
                              KTR_PARAM_MIP_OUTINTERVAL = 1,
                              KTR_PARAM_MIP_MAXNODES = 10000,
                              KTR_PARAM_HESSIAN_NO_F = KTR_HESSIAN_NO_F_ALLOW))
x_U = [2,2,1]
@variable(m, x_U[i] >= x[i=1:3] >= 0)
@variable(m, y[4:6], Bin)

@NLobjective(m, Min, 10 + 10*x[1] - 7*x[3] + 5*y[4] + 6*y[5] + 8*y[6] - 18*log(x[2]+1) - 19.2*log(x[1]-x[2]+1))
@NLconstraints(m, begin
    0.8*log(x[2] + 1) + 0.96*log(x[1] - x[2] + 1) - 0.8*x[3] >= 0
    log(x[2] + 1) + 1.2*log(x[1] - x[2] + 1) - x[3] - 2*y[6] >= -2
    x[2] - x[1] <= 0
    x[2] - 2*y[4] <= 0
    x[1] - x[2] - 2*y[5] <= 0
    y[4] + y[5] <= 1
end)
solve(m)
