from __future__ import division
import sys, os, numpy, itertools
from numpy.linalg import inv
from cvxopt import matrix, solvers
from cvxopt.modeling import op, dot, variable
from scipy.optimize import minimize

# consider two suppliers, two types
# ------------------------------------------------------------------------------
# Input
# ------------------------------------------------------------------------------
'''
So far I have tested with 1 supplier and multiple types; and with two suppliers
with a single type. Both work. I'm trying to figure out why with 2 suppliers and
2 types is not working
'''
# types and distribution
# types = {0:0.1, 1:0.2}
# fm = {0:[0.5,0.5]} # marginals of distribution
types = {0:0.1, 1:0.2}
ntypes = len(types)
fm = {0:[0.6,0.4],1:[0.75,0.25]} # both firms have only one type
# types = {0:0.1, 1:0.2}
# fm = {0:[0.5,0.5],1:[0.5,0.5]} # marginals of distribution

# suppliers
nsupp = len(fm)

# type space = cross product of types among suppliers
Theta = [p for p in itertools.product(range(ntypes), repeat=nsupp)]

# joint distribution based on marginals
f = {}
for perm in Theta:
    f[perm] = numpy.prod([fm[i][perm[i]] for i in fm])

# alphas
a_1, a_2 = 0.25, 0.5
# matrix for demand
r_1, r_2 = 1,1.1
gamma = 0.5

if r_1 + r_2 < 2*gamma:
    print '***ERROR: does not satisfy consistency check'
    sys.exit(1)

# build input matrices
if nsupp == 2:
    nalpha = numpy.array([a_1, a_2])
    nGamma = numpy.array([[r_1,-gamma],[-gamma,r_2]])
if nsupp == 1:
    nalpha = numpy.array([a_1])
    nGamma = numpy.array([[r_1]])

nD = inv(nGamma)
nc = numpy.dot(nD, nalpha)

# ------------------------------------------------------------------------------
# CENTRALIZED WITH IC AND IR
# ------------------------------------------------------------------------------
nvars = 2*len(Theta)*nsupp # here we include transfers and allocations
x0 = numpy.ones(nvars)

# this dict tell us the probability of all other subjects being of their type $f_{-i}(\theta_{-i})$
f_woi = {i:{j:0 for j in Theta} for i in range(nsupp)}
for i in f_woi:
    print i
    supp_woi = range(nsupp)
    del supp_woi[i]
    print supp_woi
    for perm in Theta:
        print perm
        f_woi[i][perm] = numpy.prod([fm[j][perm[i]] for j in supp_woi])
    print f_woi[i]

if nsupp == 1:
    f_woi = {i:{j:1 for j in Theta} for i in range(nsupp)}

# individual rationality constraints
bIR = numpy.zeros((ntypes*nsupp, nvars))
row = 0
for i in range(nsupp):
    for j in range(ntypes):
        # we assume that the type of employee i is j
        # now we find all scenarios where employee i is of type j
        idxs = [k for k in range(len(Theta)) if Theta[k][i] == j]
        print i, j, idxs
        for k in idxs:
            bIR[row, nsupp*k+i] = -types[j]*f_woi[i][Theta[k]]
            bIR[row, nsupp*(k+len(Theta))+i] = f_woi[i][Theta[k]]
        row+=1
print bIR

# incentive compatibility constraints
bIC = numpy.zeros((ntypes*(ntypes-1)*nsupp,nvars))
row = 0
for i in range(nsupp):
    for j in range(ntypes):
        other_types = range(ntypes)
        del other_types[j]
        for k in other_types:
            #row = i*ntypes*(ntypes-1) + j*(ntypes-1) + k
            print j, other_types, row
            idxs_j = [l for l in range(len(Theta)) if Theta[l][i] == j]
            idxs_k = [l for l in range(len(Theta)) if Theta[l][i] == k]
            for l in idxs_j:
                bIC[row, nsupp*l+i] = -types[j]*f_woi[i][Theta[l]]
                bIC[row, nsupp*(l+len(Theta))+i] = f_woi[i][Theta[l]]
            for l in idxs_k:
                bIC[row, nsupp*l+i] = types[j]*f_woi[i][Theta[l]]
                bIC[row, nsupp*(l+len(Theta))+i] = -f_woi[i][Theta[l]]
            row+=1
print bIC

# # non-negativity constraints for allocation
# bNN = numpy.zeros((len(Theta)*nsupp, 2*len(Theta)*nsupp)) # non-negativity constraints
# bNN[0:len(Theta)*nsupp, 0:len(Theta)*nsupp] = numpy.identity(len(Theta)*nsupp)
# print bNN

# non-negativity constraints for allocation and transfers
bNN = numpy.identity(nvars) # non-negativity constraints
print bNN


# put all inequality constraints on the same matrix
bG = bNN
bG = numpy.append(bG, bIR, axis=0)
if bIC.shape[0]>0: # if there is only one type there are no IC constraints
    bG = numpy.append(bG, bIC, axis=0)
print bG


bD = numpy.zeros((nvars, nvars))
bA = numpy.zeros((len(Theta), nvars))
bh = numpy.zeros(bG.shape[0])
bb = numpy.ones(len(Theta))
bq = numpy.zeros(nvars)

bq = numpy.zeros(nvars)
for i in range(len(Theta)):
    bD[nsupp*i:nsupp*(i+1),nsupp*i:nsupp*(i+1)] = f[Theta[i]]*nD
    bA[i,nsupp*i:nsupp*(i+1)] = numpy.ones(nsupp)
    bq[nsupp*i:nsupp*(i+1)] = f[Theta[i]]*nc
    bq[(len(Theta)*nsupp + nsupp*i):(len(Theta)*nsupp + nsupp*(i+1))] = -f[Theta[i]]*numpy.ones(nsupp)




from numpy.linalg import matrix_rank

print bA
print matrix_rank(bA)
print bG
print matrix_rank(bG)

test_matrix = numpy.append(bA, bG, axis=0)
print test_matrix
print matrix_rank(test_matrix)
# sys.exit(1)


# shapes
print 'x', x0.shape
print 'D', bD.shape
print 'q', bq.shape
print 'G', bG.shape
print 'h', bh.shape
print 'A', bA.shape
print 'b', bb.shape

# obj
D = matrix(bD)
q = matrix(bq)
# ineq constraints
G = matrix(bG, (bG.shape[0],bG.shape[1]))
h = matrix(bh)
# eq constraints
A = matrix(bA,(bA.shape[0],bA.shape[1]))
b = matrix(bb)


sol=solvers.qp(D, -q, -G, h, A, b)
print(sol['x'])


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

def EasyProblemScipy2():
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

def EasyProblemCVXOPT2():
    # now we try to implement using diagonal matrices depending on the number of scenarios
    x0 = numpy.ones(len(Theta)*nsupp)
    bD = numpy.zeros((len(x0),len(x0)))
    bA = numpy.zeros((len(Theta), len(x0)))
    bG = -numpy.identity(len(x0))
    bh = numpy.zeros(len(x0))
    bb = numpy.ones(len(Theta))
    bp = numpy.zeros(len(x0))
    bc = numpy.zeros(len(x0))

    print bD

    print Theta
    print f
    # now we build the blocks of the matrix
    for i in range(len(Theta)):
        bD[nsupp*i:nsupp*(i+1),nsupp*i:nsupp*(i+1)] = f[Theta[i]]*nD
        bA[i,nsupp*i:nsupp*(i+1)] = numpy.ones(nsupp)
        bc[nsupp*i:nsupp*(i+1)] = nc
        bp[nsupp*i:nsupp*(i+1)] = np
        print f[Theta[i]]*nD
        print f[Theta[i]]

    print bD
    print bA
    print bG
    print bh
    print bb
    print bp
    print bc



    D = matrix(bD)
    G = matrix(bG)
    h = matrix(bh)
    A = matrix(bA)
    b = matrix(bb)
    p = matrix(bp)
    c = matrix(bc)


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
