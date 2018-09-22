using JuMP
using NLopt

m = Model(solver=NLoptSolver(algorithm=:LD_MMA))

a1 = 2
b1 = 0
a2 = -1
b2 = 1

@variable(m, x1, start = 0)
@variable(m, x2 >= 0)
@variable(m, x3)

@NLobjective(m, Min, sqrt(x2))
@NLconstraint(m, x2 >= (a1*x1+b1)^3)
@NLconstraint(m, x2 >= (a2*x1+b2)^3)
# @constraint(m, x2 == x3)
@constraint(m, x2 >= x3)
@constraint(m, x2 <= x3)

setvalue(x1, 1.234)
setvalue(x2, 5.678)

status = solve(m)

println("got ", getobjectivevalue(m), " at ", [getvalue(x1),getvalue(x2)])
println("Objective value: ", getobjectivevalue(m))
println("Allocations: ", getvalue(x1))
println("Transfers: ", getvalue(x2))
println("Transfers: ", getvalue(x3))
