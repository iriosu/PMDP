# ----------------------------------
# Utilities
# ----------------------------------
function combine(a, r, n, out)
    # recursively generates all permutations to generate types
    if length(a) == r
        push!(out,tuple(deepcopy(a)...))
    else
        for i in 1:n
            push!(a,i)
            combine(a, r, n, out)
            pop!(a)
        end
    end
end

function combwithrep(r,n)
    # returns list with all types; r = number suppliers, n = number types
   out = []
   combine([], r, n, out)
   return out
end

function ComputeCumulativeDistribution(fm)
    # starting from marginal distributions, returns the cumulative distribution
    Fm = Dict(i=>[] for i in keys(fm))
    for key in keys(fm)
        for i in 1:length(fm[key])
            push!(Fm[key],sum([fm[key][j] for j in 1:i]))
        end
    end
    return Fm
end


# ----------------------------------
# Method to compute virtual costs
# ----------------------------------
function virtual_cost(theta, fm, Fm, idx_s, idx_t)
    theta_i = theta[idx_t]
    f_i, F_i = fm[idx_s], Fm[idx_s]
    if idx_t == 1
        rho_i = 1
    else
        rho_i = idx_t - 1
    end
    return theta_i + (F_i[rho_i]/f_i[idx_t])*(theta_i - theta[rho_i])
end

function interim_allocation_and_transfer(fm, Theta, nsupp, ntypes, x, t)
    # receives vector of transfers and allocations and computes the interim counterpart
    # this dict tell us the probability of all other subjects being of their type $f_{-i}(\theta_{-i})$
    f_woi = Dict(i=>Dict(j=>0.0 for j in Theta) for i in 1:nsupp)

    for i in keys(f_woi)
        supp_woi = [j for j in 1:nsupp if j!=i]
        for perm in Theta
            f_woi[i][perm] = prod([fm[j][perm[j]] for j in supp_woi])
        end
    end
    Xout = Dict(i=>Dict(j=>0.0 for j=1:ntypes) for i in 1:nsupp)
    Tout = Dict(i=>Dict(j=>0.0 for j=1:ntypes) for i in 1:nsupp)
    # individual rationality constraints
    for i in 1:nsupp
        for j in 1:ntypes
            # we assume that the type of employee i is j
            # now we find all scenarios where employee i is of type j
            idxs = [k for k in 1:length(Theta) if Theta[k][i] == j]
            for k in idxs
                Xout[i][j] = Xout[i][j] + x[nsupp*(k-1)+i]*f_woi[i][Theta[k]]
                Tout[i][j] = Tout[i][j] + t[nsupp*(k-1)+i]*f_woi[i][Theta[k]]
            end
        end
    end
    return Xout, Tout
end

# ----------------------------------
# Methods to compute analytical demands
# ----------------------------------
function demand_LM(p, a, r, gamma, idx)
    # affine demand function for 2 suppliers
    r_i, a_i, p_i = 0, 0, 0
    r_j, a_j, p_j = 0, 0, 0
    if idx == 1
        r_i, a_i, p_i = r[1], a[1], p[1]
        r_j, a_j, p_j = r[2], a[2], p[2]
    else
        r_i, a_i, p_i = r[2], a[2], p[2]
        r_j, a_j, p_j = r[1], a[1], p[1]
    end
    numerator = (r_j-gamma)*a_i - (r_i-gamma)*a_j + r_i - gamma - (r_i*r_j-gamma^2)*(p_i-p_j)
    denominator = r_i + r_j - 2*gamma
    return max(0, min( numerator/denominator, 1 )  )
end

function asssortment_HM(p, loc, delta)
    out = []
    for i in 1:length(p)
        # if p[i] < min(p[j] + delta*abs(loc[j]-loc[i]) for j=1:length(p) if j!=i)
        aux = min([p[j] + delta*abs(loc[j]-loc[i]) for j=1:length(p) if j!=i]...)
        if p[i] <= aux
            push!(out,i)
        end
    end
    return out
end

function demand_HM(p, loc, delta, idx, Q)
    # affine demand function for 2 suppliers
    if idx in Q
        if length(Q) == 1
            return 1
        end
        l_j, p_j = loc[idx], p[idx]
        if idx == min(Q...)
            idx_next = min([j for j in Q if j > idx]...)
            l_k, p_k = loc[idx_next], p[idx_next]
            aux_next = (p_k-p_j+delta*abs(l_k-l_j))/(2*delta)
            return l_j + aux_next
        elseif idx == max(Q...)
            idx_prev = max([j for j in Q if j < idx]...)
            l_i, p_i = loc[idx_prev], p[idx_prev]
            aux_prev = (p_i-p_j+delta*abs(l_j-l_i))/(2*delta)
            return (1-l_j) + aux_prev
        else
            idx_next = min([j for j in Q if j > idx]...)
            idx_prev = max([j for j in Q if j < idx]...)
            l_i, p_i = loc[idx_prev], p[idx_prev]
            l_k, p_k = loc[idx_next], p[idx_next]
            aux_next = (p_k-p_j+delta*abs(l_k-l_j))/(2*delta)
            aux_prev = (p_i-p_j+delta*abs(l_j-l_i))/(2*delta)
            return aux_prev + aux_next
        end
    else
        return 0
    end
end


# ----------------------------------
# Methods to check assumptions and conditions
# ----------------------------------
function check_vc_increasing(fm)
    # function that checks if the virtual costs are increasing. If they are, it returns them. If not, returns false
    theta = [types[i] for i in keys(types)]
    sort!(theta)
    ntypes = length(types)
    nsupp = length(fm)

    Fm = ComputeCumulativeDistribution(fm)

    v = zeros(nsupp, ntypes)

    for i in 1:nsupp
        for j in 1:ntypes
            v[i,j] = virtual_cost(theta, fm, Fm, i, j)
        end
    end

    for j in 1:(ntypes-1)
        boos = v[:,j].<v[:,j+1]
        if prod(boos) == false
            println(v)
            println("***ERROR: virtual costs are not strictly increasing")
            exit()
        end
    end
    return v
end

function check_conditions_HM(Theta, fm, loc, delta)
    # check conditions for Theorem 4.1
    V = check_vc_increasing(fm) # V[i,j] = virtual cost of supplier i for type j
    nsupp, ntypes = size(V)[1], size(V)[2]

    c_star = min([abs(loc[j+1]-loc[j]) for j=1:(nsupp-1)]...)

    # codnition 1
    boo_1 = false
    for i in 1:length(Theta)
        prob = prod([fm[j][Theta[i][j]] for j=1:nsupp])
        if prob == 0
            continue
        end
        v_c = [V[j,Theta[i][j]] for j=1:nsupp]

        boos = [abs(v_c[j+1]-v_c[j])<= 0.25*delta*abs(loc[j+1]-loc[j]) for j=1:(nsupp-1)]
        if prod(boos)
            boo_1 = true
            break
        end
    end

    # condition 2
    boo_2 = true
    for i=1:nsupp
        ntypes_i = sum([fm[i][j]>0 for j=1:ntypes])
        if ntypes_i < 3
            boo_2 = false
        end
    end

    if boo_2
        for j=1:(ntypes-1)
            boos = (V[:,j+1]-V[:,j].<=delta*c_star/8)
            if !prod(boos)
                boo_2 = false
                break
            end
        end
    end

    if !prod([boo_1, boo_2])
        println("***WARNING: the conditions on the Hotelling model to obtain same results in decentralized problem are not satisfied.")
    end
    return [boo_1, boo_2]
end

function check_conditions_LM(types, fm, a, r, gamma)
    # check conditions for Theorem 5.2
    ntypes = length(types)
    nsupp = length(fm)
    Theta = combwithrep(nsupp, ntypes)

    V = check_vc_increasing(fm) # V[i,j] = virtual cost of supplier i for type j
    nsupp, ntypes = size(V)[1], size(V)[2]

    boo_1, tsd, tnsd = false, [], [] # tsd = theta splitted demand; tnsd = theta not splitted demand
    for i=1:length(Theta)
        p_s = [V[j,Theta[i][j]] for j=1:nsupp]
        d_s = []
        for j=1:nsupp
            daux = demand_LM(p_s, a, r, gamma, j)
            push!(d_s, daux)
        end
        if prod(d_s .> 0)
            boo_1 = true
            push!(tsd, i)
        else
            push!(tnsd, i)
        end
    end

    d_star = 1e12
    if boo_1 && length(tnsd) > 0
        for i=1:length(tsd)
            vals = [abs(V[j,Theta[tsd[i]][j]]-V[j,Theta[tnsd[k]][j]]) for k=1:length(tnsd) for j=1:nsupp ]
            mx = max(vals...)
            if mx < d_star
                d_star = mx
            end
        end
    end
    # condition 2
    boo_2 = true
    for i=1:nsupp
        ntypes_i = sum([fm[i][j]>0 for j=1:ntypes])
        if ntypes_i < 3
            boo_2 = false
            break
        end
    end

    if boo_2
        for j=1:(ntypes-1)
            boos = (V[:,j+1]-V[:,j].<=d_star/2)
            if !prod(boos)
                boo_2 = false
                break
            end
        end
    end
    if !prod([boo_1, boo_2])
        println("***WARNING: the conditions on the general affine model to obtain same results in decentralized problem are not satisfied.")
    end
    return [boo_1, boo_2]
end

function check_diagonal_dominant(A)
    boos = 2*diag(abs.(A)).>sum(abs.(A),2)
    if prod(boos) == false
        println("***WARNING: the matrix provided is not strictly diagonal dominant")
    end
    return prod(boos)
end

function check_consistency_demand_matrix(A)
    if sum(diag(A)) < 2*A[1,2]
        println("***ERROR: does not satisfy consistency check")
        return false
    else
        return true
    end
end

function InputsConstraintsCentralized(types, fm)
    ntypes = length(types)
    nsupp = length(fm)
    Theta = combwithrep(nsupp, ntypes)

    # joint distribution based on marginals
    f = Dict()
    for perm in Theta
        f[perm] = prod([fm[i][perm[i]] for i=1:nsupp])
    end

    # ------------------------------------------------------------------------------
    # CENTRALIZED WITH IC AND IR
    # ------------------------------------------------------------------------------
    nvars = length(Theta)*nsupp # here we include transfers and allocations
    sts = length(Theta) # size of type space

    # this dict tell us the probability of all other subjects being of their type $f_{-i}(\theta_{-i})$
    f_woi = Dict(i=>Dict(j=>0.0 for j in Theta) for i in 1:nsupp)

    for i in keys(f_woi)
        supp_woi = [j for j in 1:nsupp if j!=i]
        for perm in Theta
            f_woi[i][perm] = prod([fm[j][perm[j]] for j in supp_woi])
        end
    end

    if nsupp == 1
        f_woi = Dict(i=>Dict(j=>1.0 for j in Theta) for i in 1:nsupp)
    end

    # individual rationality constraints
    bIR_x = zeros((ntypes*nsupp, nvars))
    bIR_t = zeros((ntypes*nsupp, nvars))
    row = 1
    for i in 1:nsupp
        for j in 1:ntypes
            # we assume that the type of employee i is j
            # now we find all scenarios where employee i is of type j
            idxs = [k for k in 1:length(Theta) if Theta[k][i] == j]
            for k in idxs
                bIR_x[row, nsupp*(k-1)+i] = -types[j]*f_woi[i][Theta[k]]
                bIR_t[row, nsupp*(k-1)+i] = f_woi[i][Theta[k]]
            end
            row=row+1
        end
    end

    bIC_x = zeros((ntypes*(ntypes-1)*nsupp,nvars))
    bIC_t = zeros((ntypes*(ntypes-1)*nsupp,nvars))
    row = 1
    for i in 1:nsupp
        for j in 1:ntypes
            other_types = [k for k in 1:ntypes if k!=j]
            for k in other_types
                # print j, other_types, row
                idxs_j = [l for l in 1:length(Theta) if Theta[l][i] == j]
                idxs_k = [l for l in 1:length(Theta) if Theta[l][i] == k]
                for l in idxs_j
                    bIC_x[row, nsupp*(l-1)+i] = -types[j]*f_woi[i][Theta[l]]
                    bIC_t[row, nsupp*(l-1)+i] = f_woi[i][Theta[l]]
                end
                for l in idxs_k
                    bIC_x[row, nsupp*(l-1)+i] = types[j]*f_woi[i][Theta[l]]
                    bIC_t[row, nsupp*(l-1)+i] = -f_woi[i][Theta[l]]
                end
                row+=1
            end
        end
    end

    bG_x = vcat(bIR_x, bIC_x)
    bG_t = vcat(bIR_t, bIC_t)

    bA = zeros((length(Theta), nvars))
    bh = zeros(size(bG_x)[1])
    bb = ones(length(Theta))

    wq_t = zeros(nvars)
    q_t = zeros(nvars)

    for i in 1:length(Theta)
        bA[i,nsupp*(i-1)+1:nsupp*i] = ones(nsupp)
        wq_t[nsupp*(i-1)+1:nsupp*i] = f[Theta[i]]*ones(nsupp)
        q_t[nsupp*(i-1)+1:nsupp*i] = ones(nsupp)
    end

    return nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, f, Theta
end

function InputsExPostConstraintsCentralized(types, fm)
    ntypes = length(types)
    nsupp = length(fm)
    Theta = combwithrep(nsupp, ntypes)

    # joint distribution based on marginals
    f = Dict()
    for perm in Theta
        f[perm] = prod([fm[i][perm[i]] for i=1:nsupp])
    end

    # ------------------------------------------------------------------------------
    # CENTRALIZED WITH IC AND IR
    # ------------------------------------------------------------------------------
    nvars = length(Theta)*nsupp # here we include transfers and allocations
    sts = length(Theta) # size of type space

    # this dict tell us the probability of all other subjects being of their type $f_{-i}(\theta_{-i})$
    f_woi = Dict(i=>Dict(j=>0.0 for j in Theta) for i in 1:nsupp)

    for i in keys(f_woi)
        supp_woi = [j for j in 1:nsupp if j!=i]
        for perm in Theta
            f_woi[i][perm] = prod([fm[j][perm[j]] for j in supp_woi])
        end
    end

    if nsupp == 1
        f_woi = Dict(i=>Dict(j=>1.0 for j in Theta) for i in 1:nsupp)
    end

    # individual rationality constraints
    bIR_x = zeros((ntypes*nsupp, nvars))
    bIR_t = zeros((ntypes*nsupp, nvars))
    row = 1
    for i in 1:nsupp
        for j in 1:ntypes
            # we assume that the type of employee i is j
            # now we find all scenarios where employee i is of type j
            idxs = [k for k in 1:length(Theta) if Theta[k][i] == j]
            for k in idxs
                bIR_x[row, nsupp*(k-1)+i] = -types[j]*f_woi[i][Theta[k]]
                bIR_t[row, nsupp*(k-1)+i] = f_woi[i][Theta[k]]
            end
            row=row+1
        end
    end

    bIC_x = zeros((ntypes*(ntypes-1)*nsupp,nvars))
    bIC_t = zeros((ntypes*(ntypes-1)*nsupp,nvars))
    row = 1
    for i in 1:nsupp
        for j in 1:ntypes
            other_types = [k for k in 1:ntypes if k!=j]
            for k in other_types
                # print j, other_types, row
                idxs_j = [l for l in 1:length(Theta) if Theta[l][i] == j]
                idxs_k = [l for l in 1:length(Theta) if Theta[l][i] == k]
                for l in idxs_j
                    bIC_x[row, nsupp*(l-1)+i] = -types[j]*f_woi[i][Theta[l]]
                    bIC_t[row, nsupp*(l-1)+i] = f_woi[i][Theta[l]]
                end
                for l in idxs_k
                    bIC_x[row, nsupp*(l-1)+i] = types[j]*f_woi[i][Theta[l]]
                    bIC_t[row, nsupp*(l-1)+i] = -f_woi[i][Theta[l]]
                end
                row+=1
            end
        end
    end

    bG_x = vcat(bIR_x, bIC_x)
    bG_t = vcat(bIR_t, bIC_t)

    bA = zeros((length(Theta), nvars))
    bh = zeros(size(bG_x)[1])
    bb = ones(length(Theta))

    wq_t = zeros(nvars)
    q_t = zeros(nvars)

    for i in 1:length(Theta)
        bA[i,nsupp*(i-1)+1:nsupp*i] = ones(nsupp)
        wq_t[nsupp*(i-1)+1:nsupp*i] = f[Theta[i]]*ones(nsupp)
        q_t[nsupp*(i-1)+1:nsupp*i] = ones(nsupp)
    end

    return nsupp, ntypes, nvars, sts, bG_x, bG_t, bh, bA, bb, wq_t, f, Theta
end
