using JuMP
using Ipopt


m = Model(solver="Ipopt")
@variable(m, x)
@variable(m, y)

setValue(x, 0.0); setValue(y, 0.0)
@NLObjective(m, Min, (1-x)^2 + 100(y-x^2)^2)

solve(m)
println("x = ", getValue(x), " y = ", getValue(y))

# adding a (linear) constraint
@addConstraint(m, x + y == 10)
solve(m)
println("x = ", getValue(x), " y = ", getValue(y))




exit()

# objective
function eval_f(x)
  return x[1]^2 + x[2]^2
end

# constraints
function eval_g(x, g)
  # Bad: g    = zeros(2)  # Allocates new array
  # OK:  g[:] = zeros(2)  # Modifies 'in place'
  g[1] = x[1]   * x[2]
end

# gradient of objective
function eval_grad_f(x, grad_f)
  # Bad: grad_f    = zeros(4)  # Allocates new array
  # OK:  grad_f[:] = zeros(4)  # Modifies 'in place'
  grad_f[1] = 2*x[1]
  grad_f[2] = 2*x[2]
end

# jacobian matrix of constraints
function eval_jac_g(x, mode, rows, cols, values)
  if mode == :Structure
    # Constraint (row) 1
    rows[1] = 1; cols[1] = 1
    rows[2] = 1; cols[2] = 2
  else
    # Constraint (row) 1
    values[1] = x[2]  # 1,1
    values[2] = x[1]  # 1,2
  end
end

# hessian matrix  is the sum of the hessian of objective plus the constraints times multipliers
function eval_h(x, mode, rows, cols, obj_factor, lambda, values)
  if mode == :Structure
    # Symmetric matrix, fill the lower left triangle only
    idx = 1
    for row = 1:2
      for col = 1:row
        rows[idx] = row
        cols[idx] = col
        idx += 1
      end
    end
  else
    # Again, only lower left triangle
    # Objective
    values[1] = obj_factor * 2 # 1,1
    values[3] = obj_factor * 2 # 2,2

    # First constraint
    values[2] += lambda[1]   # 2,1
  end
end

n = 2
x_L = [0.0,0.0]
x_U = [5.0, 10.0]

m = 1
g_L = [1.0]
g_U = [1.0]

prob = createProblem(n, x_L, x_U, m, g_L, g_U, 2, 3,
                     eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h)

prob.x = [0.1, 10]
status = solveProblem(prob)

println(Ipopt.ApplicationReturnStatus[status])
println(prob.x)
println(prob.obj_val)


exit()












# hs071
# min x1 * x4 * (x1 + x2 + x3) + x3
# st  x1 * x2 * x3 * x4 >= 25
#     x1^2 + x2^2 + x3^2 + x4^2 = 40
#     1 <= x1, x2, x3, x4 <= 5
# Start at (1,5,5,1)
# End at (1.000..., 4.743..., 3.821..., 1.379...)

# objective
function eval_f(x)
  return x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3]
end

# constraints
function eval_g(x, g)
  # Bad: g    = zeros(2)  # Allocates new array
  # OK:  g[:] = zeros(2)  # Modifies 'in place'
  g[1] = x[1]   * x[2]   * x[3]   * x[4]
  g[2] = x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2
end

# gradient of objective
function eval_grad_f(x, grad_f)
  # Bad: grad_f    = zeros(4)  # Allocates new array
  # OK:  grad_f[:] = zeros(4)  # Modifies 'in place'
  grad_f[1] = x[1] * x[4] + x[4] * (x[1] + x[2] + x[3])
  grad_f[2] = x[1] * x[4]
  grad_f[3] = x[1] * x[4] + 1
  grad_f[4] = x[1] * (x[1] + x[2] + x[3])
end

# jacobian matrix of constraints
function eval_jac_g(x, mode, rows, cols, values)
  if mode == :Structure
    # Constraint (row) 1
    rows[1] = 1; cols[1] = 1
    rows[2] = 1; cols[2] = 2
    rows[3] = 1; cols[3] = 3
    rows[4] = 1; cols[4] = 4
    # Constraint (row) 2
    rows[5] = 2; cols[5] = 1
    rows[6] = 2; cols[6] = 2
    rows[7] = 2; cols[7] = 3
    rows[8] = 2; cols[8] = 4
  else
    # Constraint (row) 1
    values[1] = x[2]*x[3]*x[4]  # 1,1
    values[2] = x[1]*x[3]*x[4]  # 1,2
    values[3] = x[1]*x[2]*x[4]  # 1,3
    values[4] = x[1]*x[2]*x[3]  # 1,4
    # Constraint (row) 2
    values[5] = 2*x[1]  # 2,1
    values[6] = 2*x[2]  # 2,2
    values[7] = 2*x[3]  # 2,3
    values[8] = 2*x[4]  # 2,4
  end
end

# hessian matrix  is the sum of the hessian of objective plus the constraints times multipliers
function eval_h(x, mode, rows, cols, obj_factor, lambda, values)
  if mode == :Structure
    # Symmetric matrix, fill the lower left triangle only
    idx = 1
    for row = 1:4
      for col = 1:row
        rows[idx] = row
        cols[idx] = col
        idx += 1
      end
    end
  else
    # Again, only lower left triangle
    # Objective
    values[1] = obj_factor * (2*x[4])  # 1,1
    values[2] = obj_factor * (  x[4])  # 2,1
    values[3] = 0                      # 2,2
    values[4] = obj_factor * (  x[4])  # 3,1
    values[5] = 0                      # 3,2
    values[6] = 0                      # 3,3
    values[7] = obj_factor * (2*x[1] + x[2] + x[3])  # 4,1
    values[8] = obj_factor * (  x[1])  # 4,2
    values[9] = obj_factor * (  x[1])  # 4,3
    values[10] = 0                     # 4,4

    # First constraint
    values[2] += lambda[1] * (x[3] * x[4])  # 2,1
    values[4] += lambda[1] * (x[2] * x[4])  # 3,1
    values[5] += lambda[1] * (x[1] * x[4])  # 3,2
    values[7] += lambda[1] * (x[2] * x[3])  # 4,1
    values[8] += lambda[1] * (x[1] * x[3])  # 4,2
    values[9] += lambda[1] * (x[1] * x[2])  # 4,3

    # Second constraint
    values[1]  += lambda[2] * 2  # 1,1
    values[3]  += lambda[2] * 2  # 2,2
    values[6]  += lambda[2] * 2  # 3,3
    values[10] += lambda[2] * 2  # 4,4
  end
end

n = 4
x_L = [1.0, 1.0, 1.0, 1.0]
x_U = [5.0, 5.0, 5.0, 5.0]

m = 2
g_L = [25.0, 40.0]
g_U = [2.0e19, 40.0]

###
# createProblem
#
# (C function: CreateIpoptProblem):
#
# function createProblem(
#   n::Int,                     # Number of variables
#   x_L::Vector{Float64},       # Variable lower bounds
#   x_U::Vector{Float64},       # Variable upper bounds
#   m::Int,                     # Number of constraints that are not bounds
#   g_L::Vector{Float64},       # Constraint lower bounds
#   g_U::Vector{Float64},       # Constraint upper bounds ==> equality constrained defined as lower=upper bound
#   nele_jac::Int,              # Number of non-zeros in Jacobian
#   nele_hess::Int,             # Number of non-zeros in Hessian
#   eval_f,                     # Callback: objective function
#   eval_g,                     # Callback: constraint evaluation
#   eval_grad_f,                # Callback: objective function gradient
#   eval_jac_g,                 # Callback: Jacobian evaluation
#   eval_h = nothing)           # Callback: Hessian evaluation
###

prob = createProblem(n, x_L, x_U, m, g_L, g_U, 8, 10,
                     eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h)

prob.x = [1.0, 5.0, 5.0, 1.0]
status = solveProblem(prob)

println(Ipopt.ApplicationReturnStatus[status])
println(prob.x)
println(prob.obj_val)
