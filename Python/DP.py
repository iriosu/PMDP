from __future__ import division
import sys, os, numpy, itertools, math
from numpy.linalg import inv
from cvxopt import matrix, solvers, spdiag
from cvxopt.modeling import op, dot, variable
from scipy.optimize import minimize



'''
Solves the problem of Decentralized Procurement
'''
def GenerateAllQ():
    def findsubsets(S,m):
        return set(itertools.combinations(S, m))

    suppliers = range(2)
    Q=[]
    for i in range(1,len(suppliers)+1):
        aux = findsubsets(suppliers,i)
        for j in aux:
            Q.append(j)
            print i, j

    print Q

    M = [p for p in itertools.product(range(3), repeat=4)]
    print M
    print len(M)

    oQ = []
    for i in range(len(M)):
        ra = []
        for j in range(len(M[i])):
            if M[i][j] == 0:
                ra.extend([1,0])
            elif M[i][j] == 1:
                ra.extend([0,1])
            else:
                ra.extend([1,1])
        ra.extend(ra)
        oQ.append(ra)

    print oQ
    print len(oQ)





def GenerateInputsOPT(types, fm, a_1, a_2, r_1, r_2, gamma):
    # suppliers
    nsupp = len(fm)
    ntypes = len(types)
    # type space = cross product of types among suppliers
    Theta = [p for p in itertools.product(range(ntypes), repeat=nsupp)]

    # joint distribution based on marginals
    f = {}
    for perm in Theta:
        f[perm] = numpy.prod([fm[i][perm[i]] for i in fm])

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
    nvars = 4*len(Theta)*nsupp # here we include allocations, prices, and dual variables
    x0 = matrix(numpy.ones(nvars))

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
    bNN = numpy.zeros((nvars,nvars)) # non-negativity constraints
    bNN[0:3*len(Theta)*nsupp, 0:3*len(Theta)*nsupp] = numpy.identity(3*len(Theta)*nsupp)

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

    print bA


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

    return D,q,G,h,A,b,nvars


def F(x = None, z = None, verbose = False):
    if x is None:  return 0, matrix(0.0, (nvars,1))

    nss = int(nvars/4)
    f = numpy.zeros(nss+1)
    obj = 0.5 * x0.T * D * x0 + q.T*x0
    f[0] = obj[0]
    for i in range(nss):
        f[i+1] = x0[i]*x0[2*nss+i]

    print f

    Df = numpy.zeros((nss+1,nvars))
    D_x = D*x0
    print D_x.T
    print q
    Df[0,:] = (D_x+q).T
    for i in range(nss):
        Df[i+1,i] = x0[2*nss+i]
        Df[i+1,2*nss+i] = x0[i]

    f = matrix(f)
    Df = matrix(Df)

    if z is None:  return f, Df
    Hs = {i:numpy.zeros((nvars, nvars)) for i in range(nss+1)}
    Hs[0] = D
    for i in range(nss):
        Hs[i+1][i,2*nss+i] = 1
        Hs[i+1][2*nss+i,i] = 1

    z = numpy.ones(nss+1)
    H = numpy.zeros((nvars, nvars))
    for i in range(nss+1):
        H = H + z[i]*Hs[i]

    print H
    H = matrix(H)
    return f, Df, H










if __name__ == '__main__':
    # consider two suppliers, two types
    # ------------------------------------------------------------------------------
    # Input
    # ------------------------------------------------------------------------------

    # 2 types, 1 supplier
    # types = {0:0.1, 1:0.2}
    # fm = {0:[0.5,0.5]} # marginals of distribution

    # 2 types, 2 suppliers
    types = {0:0.1, 1:0.2}

    fm = {0:[0.6,0.4],1:[0.75,0.25]} # both firms have only one type

    # 3 types, 2 suppliers
    # types = {0:0.1, 1:0.2, 2:0.25}
    # ntypes = len(types)
    # fm = {0:[0.6,0.2, 0.2],1:[0.5,0.25, 0.25]} # both firms have only one type

    # alphas
    a_1, a_2 = 0.25, 0.5
    # matrix for demand
    r_1, r_2 = 1,1.1
    gamma = 0.5

    D,q,G,h,A,b,nvars = GenerateInputsOPT(types, fm, a_1, a_2, r_1, r_2, gamma)
    # sol=solvers.qp(D, -q, -G, h, A, b)
    # print(sol['x'])

    x0 = matrix(numpy.ones(nvars))
    z0 = matrix(numpy.ones(int(nvars/4)+1))
    f,Df,H = F(x0,z0)
    print f
    print Df
    print H


    dims = {'l': h.size[0], 'q': [], 's': []}
    print dims
    # sol = solvers.cp(F, G, h, dims)
    # print(sol['x'])

    #
    # f = numpy.zeros(len(Theta)*nsupp+1)
    # obj = 0.5 * x0.T * D * x0 + q.T*x0
    # f[0] = obj[0]
    # for i in range(len(Theta)*nsupp):
    #     f[i+1] = x0[i]*x0[2*len(Theta)*nsupp+i]
    # print f
    #
    # Df = numpy.zeros((len(Theta)*nsupp+1,nvars))
    # D_x = D*x0
    # print D_x.T
    # print q
    # Df[0,:] = (D_x+q).T
    # for i in range(len(Theta)*nsupp):
    #     Df[i+1,i] = x0[2*len(Theta)*nsupp+i]
    #     Df[i+1,2*len(Theta)*nsupp+i] = x0[i]
    #
    # print Df
    #
    # Hs = {i:numpy.zeros((nvars, nvars)) for i in range(len(Theta)*nsupp+1)}
    # Hs[0] = bD
    # for i in range(len(Theta)*nsupp):
    #     Hs[i+1][i,2*len(Theta)*nsupp+i] = 1
    #     Hs[i+1][2*len(Theta)*nsupp+i,i] = 1
    #
    # z = numpy.ones(len(Theta)*nsupp+1)
    # H = numpy.zeros((nvars, nvars))
    # for i in range(len(Theta)*nsupp+1):
    #     H = H + z[i]*Hs[i]
    #
    # print H
    # H = matrix(H)
    #
    #
    # sys.exit(1)
