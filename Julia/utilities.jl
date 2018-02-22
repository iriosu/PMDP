function combine(a, r, n, out)
    # recursively generates all permutations to generate types
    if length(a) == r
        println(a)
        println(tuple(a))
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
