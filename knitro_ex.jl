using KNITRO, JuMP

## Solve test problem 1 (Synthesis of processing system) in
 #  M. Duran & I.E. Grossmann, "An outer approximation algorithm for
 #  a class of mixed integer nonlinear programs", Mathematical
 #  Programming 36, pp. 307-339, 1986.  The problem also appears as
 #  problem synthes1 in the MacMINLP test set.


 m = Model(solver=KnitroSolver(mip_method = KTR_MIP_METHOD_BB,
                               ms_enable = 1,
                               algorithm = KTR_ALG_ACT_CG,
                               outmode = KTR_OUTMODE_SCREEN,
                               KTR_PARAM_OUTLEV = KTR_OUTLEV_ALL,
                               KTR_PARAM_MIP_OUTINTERVAL = 1,
                               KTR_PARAM_MIP_MAXNODES = 10000,
                               KTR_PARAM_HESSIAN_NO_F = KTR_HESSIAN_NO_F_ALLOW))
 x_U = [5,10,10,5]
 @variable(m, x_U[i] >= x[i=1:4] >= 0)


 @NLobjective(m, Max, x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2)
 # @constraint(m,drel[i in 1:2],x[2*(i-1)+1]+x[2*(i-1)+2] == 1)
 @NLconstraint(m,nll[i in 1:2],x[2*(i-1)+1]*x[2*(i-1)+2] == 1)
 solve(m)




exit()


 m = Model(solver=KnitroSolver(mip_method = KTR_MIP_METHOD_BB,
                               ms_enable = 1,
                               algorithm = KTR_ALG_ACT_CG,
                               outmode = KTR_OUTMODE_SCREEN,
                               KTR_PARAM_OUTLEV = KTR_OUTLEV_ALL,
                               KTR_PARAM_MIP_OUTINTERVAL = 1,
                               KTR_PARAM_MIP_MAXNODES = 10000,
                               KTR_PARAM_HESSIAN_NO_F = KTR_HESSIAN_NO_F_ALLOW))
 x_U = [5,10]
 @variable(m, x_U[i] >= x[i=1:2] >= 0)


 @NLobjective(m, Max, x[1]^2 + x[2]^2)
 @NLconstraints(m, begin
     x[1]*x[2] >= 1
     x[1]*x[2] <= 1
 end)
 solve(m)




exit()
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
