
include("utilities.jl")
function ReadOutputHM(infile)
    f = open(infile)
    lines = readlines(f)
    out = Dict()
    for l in lines
        pieces = split(l, ";")
        delta = parse(Float64, pieces[1])
        out[delta] = Dict()
        out[delta]["obj_cent"] = round(parse(Float64, pieces[2]),4)
        out[delta]["x_cent"] = [round(parse(Float64, pieces[i]),4) for i=3:10]
        out[delta]["t_cent"] = [round(parse(Float64, pieces[i]),4) for i=11:18]

        out[delta]["obj_dec"] = round(parse(Float64, pieces[19]),4)
        out[delta]["x_dec"] = [round(parse(Float64, pieces[i]),4) for i=20:27]
        out[delta]["t_dec"] = [round(parse(Float64, pieces[i]),4) for i=28:35]
    end
    return out
end

function CheckBindingIC(infile, types)

    sol = ReadOutputHM(infile)
    substr = split(infile, "HM")
    subsubstr = split(substr[2], "_")
    distr = [parse(Float64,subsubstr[i]) for i=2:3]
    fm = Dict(1=>distr,2=>distr)
    nsupp, ntypes, nvars, sts, G_x, G_t, h, A, b, q_t, f, Theta = InputsConstraintsCentralized(types, fm)
    # G_x, G_t, h
    for k in keys(sol)
        boos = abs.(G_x*sol[k]["x_cent"] + G_t*sol[k]["t_cent"]-h) .<=1e-2
        println(k, " ", boos)
        # println(G_x*sol[k]["x_cent"] + G_t*sol[k]["t_cent"])
    end
end

function CheckExPostIR(infile, types)
    sol = ReadOutputHM(infile)
    substr = split(infile, "HM")
    subsubstr = split(substr[2], "_")
    distr = [parse(Float64,subsubstr[i]) for i=2:3]
    fm = Dict(1=>distr,2=>distr)
    nsupp, ntypes, nvars, sts, G_x, G_t, h, A, b, q_t, f, Theta = InputsConstraintsCentralized(types, fm)

    outprob = Dict(k=>Dict(i=>Dict("cent"=>0.0, "dec"=>0.0) for i=1:nsupp) for k in keys(sol))
    outprof = Dict(k=>Dict(i=>Dict("cent"=>0.0, "dec"=>0.0) for i=1:nsupp) for k in keys(sol))
    for k in keys(sol)
        println("Processing delta = ", k)
        for s=1:length(Theta)
            theta = [types[Theta[s][i]] for i=1:nsupp]
            for i=1:nsupp
                u_cent = sol[k]["t_cent"][nsupp*(s-1)+i] - theta[i]*sol[k]["x_cent"][nsupp*(s-1)+i]
                u_dec = sol[k]["t_dec"][nsupp*(s-1)+i] - theta[i]*sol[k]["x_dec"][nsupp*(s-1)+i]
                # println(u_cent)
                # println(Theta[s])
                # println(f[Theta[s]])
                if u_cent < 0
                    outprob[k][i]["cent"] += f[Theta[s]]
                    outprof[k][i]["cent"] += u_cent*f[Theta[s]]
                end
                if u_dec < 0
                    outprob[k][i]["dec"] += f[Theta[s]]
                    outprof[k][i]["dec"] += u_dec*f[Theta[s]]
                end
            end
        end
    end
    return outprob, outprof
end

function CheckExAnteIR(infile, types)
    sol = ReadOutputHM(infile)
    substr = split(infile, "HM")
    subsubstr = split(substr[2], "_")
    distr = [parse(Float64,subsubstr[i]) for i=2:3]
    fm = Dict(1=>distr,2=>distr)
    nsupp, ntypes, nvars, sts, G_x, G_t, h, A, b, q_t, f, Theta = InputsConstraintsCentralized(types, fm)

    outprof = Dict(k=>Dict(i=>Dict("cent"=>0.0, "dec"=>0.0) for i=1:nsupp) for k in keys(sol))
    for k in keys(sol)
        for s=1:length(Theta)
            theta = [types[Theta[s][i]] for i=1:nsupp]
            for i=1:nsupp
                u_cent = sol[k]["t_cent"][nsupp*(s-1)+i] - theta[i]*sol[k]["x_cent"][nsupp*(s-1)+i]
                u_dec = sol[k]["t_dec"][nsupp*(s-1)+i] - theta[i]*sol[k]["x_dec"][nsupp*(s-1)+i]
                # println(u_cent)
                # println(Theta[s])
                # println(f[Theta[s]])
                outprof[k][i]["cent"] += u_cent*f[Theta[s]]
                outprof[k][i]["dec"] += u_dec*f[Theta[s]]
            end
        end
    end
    return outprof
end

function GetConsecutivePairs(nsupp)
    out = []
    for i=1:(nsupp-1)
        push!(out, (i,i+1))
        push!(out, (i+1,i))
    end
    return out
end

function Elasticities(infile, types)
    sol = ReadOutputHM(infile)
    substr = split(infile, "HM")
    subsubstr = split(substr[2], "_")
    distr = [parse(Float64,subsubstr[i]) for i=2:3]
    fm = Dict(1=>distr,2=>distr)
    nsupp, ntypes, nvars, sts, G_x, G_t, h, A, b, q_t, f, Theta = InputsConstraintsCentralized(types, fm)
    pairs = GetConsecutivePairs(nsupp)
    own_e = Dict(k=>Dict(i=>Dict("cent"=>0.0, "dec"=>0.0) for i=1:nsupp) for k in keys(sol))
    cross_e = Dict(k=>Dict(pairs[i]=>Dict("cent"=>0.0, "dec"=>0.0) for i=1:length(pairs)) for k in keys(sol))
    own_p = Dict(k=>Dict(i=>Dict("cent"=>0.0, "dec"=>0.0) for i=1:nsupp) for k in keys(sol))
    cross_p = Dict(k=>Dict(pairs[i]=>Dict("cent"=>0.0, "dec"=>0.0) for i=1:length(pairs)) for k in keys(sol))

    for k in keys(sol)
        println("Processing delta = ", k)
        for s=1:length(Theta)
            theta = [types[Theta[s][i]] for i=1:nsupp]
            for i=1:nsupp
                if sol[k]["x_cent"][nsupp*(s-1)+i] > 1e-2
                    own_e[k][i]["cent"] += -f[Theta[s]] * 0.5 * k * sol[k]["t_cent"][nsupp*(s-1)+i]/(sol[k]["x_cent"][nsupp*(s-1)+i]^2)
                    own_p[k][i]["cent"] += f[Theta[s]]
                end
                if sol[k]["x_dec"][nsupp*(s-1)+i] > 1e-2
                    own_e[k][i]["dec"] += -f[Theta[s]] * 0.5 * k * sol[k]["t_dec"][nsupp*(s-1)+i]/(sol[k]["x_dec"][nsupp*(s-1)+i]^2)
                    own_p[k][i]["dec"] += f[Theta[s]]
                end
            end
            for r=1:length(pairs)
                r1, r2 = pairs[r][1], pairs[r][2]
                if sol[k]["x_cent"][nsupp*(s-1)+r1] > 1e-2 && sol[k]["x_cent"][nsupp*(s-1)+r2] > 1e-2
                    cross_e[k][pairs[r]]["cent"] += f[Theta[s]] * 0.5 * k * sol[k]["t_cent"][nsupp*(s-1)+r2]/(sol[k]["x_cent"][nsupp*(s-1)+r1]*sol[k]["x_cent"][nsupp*(s-1)+r2])
                    cross_p[k][pairs[r]]["cent"] += f[Theta[s]]
                end
                if sol[k]["x_dec"][nsupp*(s-1)+r1] > 1e-2 && sol[k]["x_dec"][nsupp*(s-1)+r2] > 1e-2
                    cross_e[k][pairs[r]]["dec"] += f[Theta[s]] * 0.5 * k * sol[k]["t_dec"][nsupp*(s-1)+r2]/(sol[k]["x_dec"][nsupp*(s-1)+r1]*sol[k]["x_dec"][nsupp*(s-1)+r2])
                    cross_p[k][pairs[r]]["dec"] += f[Theta[s]]
                end
            end
        end
    end
    return own_e, cross_e, own_p, cross_p
end

#### MAIN ####
types = Dict(1=>10, 2=>12)
# out = ReadOutputHM("outputs/simulations_outcome_HM_0.1_0.9_inelastic.txt")
# CheckBindingIC("outputs/simulations_outcome_HM_0.1_0.9_inelastic.txt", types)
# outprob = CheckExPostIR("outputs/simulations_outcome_HM_0.4_0.6_inelastic.txt", types)
# for k in keys(outprob)
#     println("Printing for delta = ", k)
#     for i in keys(outprob[k])
#         println("    Supplier ", i)
#         println("        Cent: ", outprob[k][i]["cent"], " Dec: ", outprob[k][i]["dec"])
#     end
# end

own_e, cross_e, own_p, cross_p = Elasticities("outputs/simulations_outcome_HM_0.4_0.6_inelastic.txt", types)
for k in keys(own_e)
    for i in keys(own_e[k])
        println(k, " ", i, " ", "cent", " ", own_e[k][i]["cent"])
    end
end

# for k in keys(out)
#     println(k, " ", out[k]["obj_cent"], " ", out[k]["obj_dec"])
#     println(out[k]["x_cent"])
#     println(out[k]["x_dec"])
#     println(out[k]["t_cent"])
#     println(out[k]["t_dec"])
# end
