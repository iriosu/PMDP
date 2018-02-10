from __future__ import division
import sys, os, numpy, itertools
from numpy.linalg import inv
from cvxopt import matrix, solvers
from cvxopt.modeling import op, dot, variable
from scipy.optimize import minimize

'''
Solves the problem of Centralized Procurement
'''
# consider two suppliers, two types
# ------------------------------------------------------------------------------
# Input
# ------------------------------------------------------------------------------

# 2 types, 1 supplier
# types = {0:0.1, 1:0.2}
# fm = {0:[0.5,0.5]} # marginals of distribution

# 2 types, 2 suppliers
types = {0:0.1, 1:0.2}
ntypes = len(types)
fm = {0:[0.6,0.4],1:[0.75,0.25]} # both firms have only one type

# 3 types, 2 suppliers
# types = {0:0.1, 1:0.2, 2:0.25}
# ntypes = len(types)
# fm = {0:[0.6,0.2, 0.2],1:[0.5,0.25, 0.25]} # both firms have only one type

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
