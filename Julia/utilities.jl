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

function virtual_cost(theta, fm, Fm, idx_s, idx_t)
    theta_i = theta[idx_t]
    f_i, F_i = fm[idx_s], Fm[idx_s]
    if idx_t == 1
        rho_i = 1
    else
        rho_i = idx_t - 1
    end
    # println("    ", theta_i, " ", F_i[rho_i], " ", f_i[idx_t], " ", theta[rho_i], " ", rho_i)
    return theta_i + (F_i[rho_i]/f_i[idx_t])*(theta_i - theta[rho_i])
end

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
    c_star = min([abs(loc[j+1]-loc[j]) for j=1:(nsupp-1)]...)
    V = check_vc_increasing(fm) # V[i,j] = virtual cost of supplier i for type j
    boo = false
    println(size(V)[2])

    for i in 1:length(Theta)
        prob = prod([fm[j][Theta[i][j]] for j=1:nsupp])
        if prob == 0
            continue
        end
        v_c = [V[j,Theta[i][j]] for j=1:nsupp]

        boos = [abs(v_c[j+1]-v_c[j])<= 0.25*delta*abs(loc[j+1]-loc[j]) for j=1:(nsupp-1)]
        println(boos)
        println(prod(boos))
    end

end
