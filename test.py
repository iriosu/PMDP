from __future__ import division
import sys, os, numpy, itertools
from numpy.linalg import inv
from cvxopt import matrix, solvers
from cvxopt.modeling import op, dot, variable
from scipy.optimize import minimize

# consider two suppliers, two types
# input types -----------
ntypes = 2
nsupp = 2
Theta = [p for p in itertools.product(range(ntypes), repeat=nsupp)]
fm = {0:[0.5,0.5],1:[0.5,0.5]} # both firms have only one type
f = {}
print Theta
for perm in Theta:
    f[perm] = numpy.prod([fm[i][perm[i]] for i in fm])
# ------------------------
a_1, a_2 = 0.5, 0.25
r_1, r_2 = 1,1
gamma = 0.5
p_1, p_2 = 0.5,0.25
if r_1 + r_2 < 2*gamma:
    print '***ERROR: does not satisfy consistency check'
    sys.exit(1)

nalpha = numpy.array([a_1, a_2])
nGamma = numpy.array([[r_1,-gamma],[-gamma,r_2]])
nD = inv(nGamma)
nc = numpy.dot(nD, nalpha)
np = numpy.array([p_1, p_2])
ne = numpy.ones(2)

def loss(x):
    fobj = 0
    ns = nD.shape[0]
    for i in range(len(Theta)):
        idi, idf = ns*i, ns*(i+1)
        fobj += f[Theta[i]]*(0.5*numpy.dot(numpy.dot(x[idi:idf], nD), x[idi:idf]) + numpy.dot((np-nc),x[idi:idf]))
    return fobj




x0 = numpy.ones(len(Theta)*nD.shape[0])
# cons = [{'type': 'eq',
#          'fun' : lambda x: numpy.sum(x[i:(i+nD.shape[0])]) - 1.0} for i in range(0,len(x0),nD.shape[0])]


cons = ({'type': 'eq','fun' : lambda x: x[0]+x[1] - 1.0},
        {'type': 'eq','fun' : lambda x: x[2]+x[3] - 1.0},
        {'type': 'eq','fun' : lambda x: x[4]+x[5] - 1.0},
        {'type': 'eq','fun' : lambda x: x[6]+x[7] - 1.0})
res = minimize(loss, x0, method='SLSQP', constraints=cons,
               bounds=[(0, numpy.inf) for i in range(len(Theta)*nD.shape[0])], options={'disp': True})

print(res.x)


sys.exit(1)
#-------------------
# easier problem
#-------------------
def EasyProblemScipy():
    def loss(x):
        return 0.5*numpy.dot(numpy.dot(x, nD), x) + numpy.dot((np-nc),x)

    cons = ({'type': 'eq',
             'fun' : lambda x: numpy.sum(x) - 1.0})

    x0 = numpy.ones(nD.shape[0])

    res = minimize(loss, x0, method='SLSQP', constraints=cons,
                   bounds=[(0, numpy.inf) for i in range(nD.shape[0])], options={'disp': True})

    print(res.x)

def EasyProblemCVXOPT():
    a_1, a_2 = 0.5, 0.25
    r_1, r_2 = 1,1
    gamma = 0.5
    p_1, p_2 = 0.5,0.25
    if r_1 + r_2 < 2*gamma:
        print '***ERROR: does not satisfy consistency check'
        sys.exit(1)

    nalpha = numpy.array([a_1, a_2])
    nGamma = numpy.array([[r_1,-gamma],[-gamma,r_2]])
    nD = inv(nGamma)
    nc = numpy.dot(nD, nalpha)
    np = numpy.array([p_1, p_2])
    ne = numpy.ones(2)

    # define cvxopt equivalents
    c = matrix(nc)
    D = matrix(nD)
    p = matrix(np)
    e = matrix(ne)

    print c-p

    G = matrix([[-1.0,0.0],[0.0,-1.0]])
    h = matrix([0.0,0.0])
    A = matrix([1.0, 1.0], (1,2))
    b = matrix(1.0)
    # note: this method minimizes; since we want to maximize, we minimize the negative of our objective
    sol=solvers.qp(D, p-c, G, h, A, b)
    print(sol['x'])




a_1, a_2 = 0.5, 0.25
r_1, r_2 = 1,1
gamma = 0.5
p_1, p_2 = 0.5,0.25
if r_1 + r_2 < 2*gamma:
    print '***ERROR: does not satisfy consistency check'
    sys.exit(1)

nalpha = numpy.array([a_1, a_2])
nGamma = numpy.array([[r_1,-gamma],[-gamma,r_2]])
nD = inv(nGamma)
nc = numpy.dot(nD, nalpha)
np = numpy.array([p_1, p_2])
ne = numpy.ones(2)

# define cvxopt equivalents
c = matrix(nc)
D = matrix(nD)
p = matrix(np)
e = matrix(ne)

G = matrix([[-1.0,0.0],[0.0,-1.0]])
h = matrix([0.0,0.0])
A = matrix([1.0, 1.0], (1,2))
b = matrix(1.0)

fobj = []
ineq = []
for perm in Theta:
    a = f[perm]
    ineq.append(( G*x[perm] <= h ))
    ineq.append(( A*x[perm] == b ))
    pc_x_a = matrix(a*(p-c),(2,1))
    fobj.append(dot(pc_x_a,x[perm]))

obj = sum(fobj)
lp2 = solvers.qp(obj, ineq)
lp2.solve()
print(lp2.objective.value())

# sol=solvers.qp(D, p-c, G, h, A, b)
# print(sol['x'])





# A = matrix([[2.,1.,-1.,0.], [1.,2.,0.,-1.]])
# b = matrix([3.,3.,0.,0.])
# c = matrix([-4.,-5.])
#
# x={}
# for perm in Theta:
#     f[perm] = numpy.prod([fm[i][perm[i]] for i in fm])
#     x[perm] = variable(2)
#
# fobj = []
# ineq = []
# for perm in Theta:
#     a = f[perm]
#     ineq.append(( A*x[perm] <= b ))
#     c_x_a = matrix(a*c,(2,1))
#     fobj.append(dot(c_x_a,x[perm]))
#
# obj = sum(fobj)
# lp2 = op(obj, ineq)
# lp2.solve()
# print(lp2.objective.value())
