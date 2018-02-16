import sys, numpy
from cvxopt import matrix, log, div, spdiag, solvers

def F(x = None, z = None, verbose = False):
    if x is None:  return 2, matrix(0.0, (2,1))
    u = x**2
    f = matrix([-sum(u), x[0]*x[1]-1, 1-x[0]*x[1]])
    Df = matrix([ [2*x[0], 2*x[1]],[x[1],x[0]],[-x[1],-x[0]] ]).T
    if z is None:  return f, Df
    print z
    H = matrix([ [2*z[0], z[1]-z[2]],[z[1]-z[2], 2*z[0]] ])
    return f, Df, H

G = matrix([ [1,0],[0,1],[-1,0],[0,-1] ], (4,2))
h = numpy.zeros(4)
h[0]=10
h[1]=5

I = numpy.identity(2)
I = numpy.append(I, -I, axis=0)
h = matrix(h)
G = matrix(I)
# print G
# print h

x0 = matrix(numpy.ones(2))
x0[1]=0.1
x0[0]=10
z0 = matrix(numpy.ones(3))

print G*x0

f,Df,H = F(x0,z0)
print f
print Df
print H

sol = solvers.cp(F, G, h)
print(sol['sl'])

# F(sol['x'], None, True)
