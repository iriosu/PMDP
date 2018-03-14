using PyPlot

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
    # to compute expected loss
    outprob = Dict(k=>Dict(i=>Dict("cent"=>0.0, "dec"=>0.0) for i=1:nsupp) for k in keys(sol))
    outprof = Dict(k=>Dict(i=>Dict("cent"=>0.0, "dec"=>0.0) for i=1:nsupp) for k in keys(sol))
    ut = Dict(k=>Dict(s=>Dict(i=>Dict("cent"=>0.0, "dec"=>0.0) for i=1:nsupp) for s=1:length(Theta)) for k in keys(sol))

    for k in keys(sol)
        println("Processing delta = ", k)
        for s=1:length(Theta)
            theta = [types[Theta[s][i]] for i=1:nsupp]
            for i=1:nsupp
                u_cent = sol[k]["t_cent"][nsupp*(s-1)+i] - theta[i]*sol[k]["x_cent"][nsupp*(s-1)+i]
                u_dec = sol[k]["t_dec"][nsupp*(s-1)+i] - theta[i]*sol[k]["x_dec"][nsupp*(s-1)+i]
                ut[k][s][i]["cent"] = u_cent
                ut[k][s][i]["dec"] = u_dec
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
    return ut, outprob, outprof
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
        for s=1:length(Theta)
            theta = [types[Theta[s][i]] for i=1:nsupp]
            for i=1:nsupp
                println(k ," ", i, " ", s, " ", sol[k]["x_cent"][nsupp*(s-1)+i])
                println(k ," ", i, " ", s, " ", sol[k]["t_cent"][nsupp*(s-1)+i])
                # println(k ," ", i, " ", s, " ", nsupp*(s-1)+i, " ", sol[k]["x_cent"])
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

function PlotObjectives(distr, elastic, version, expostir=false)
    els = nothing
    if elastic
        els = "elastic"
    else
        els = "inelastic"
    end
    for k=1:length(distr)
        dst = distr[k]
        println("========================")
        println("Distribution: ", dst)
        println("========================")
        strdist = join(dst,'_')
        if expostir == true
            infile = string("outputs/simulations_outcome_HM_", join(dst,'_'), "_expostir")
        else
            infile = string("outputs/simulations_outcome_HM_", join(dst,'_'))
        end


        if elastic
            els = "elastic"
            infile = string(infile, "_elastic.txt")
        else
            els = "inelastic"
            infile = string(infile, "_inelastic.txt")
        end
        inls = ReadOutputHM(infile)
        ds = sort([i for i in keys(inls)])
        println([inls[s]["obj_cent"] for s in ds])
        println([inls[s]["obj_dec"] for s in ds])
        if version == "centralized"
            plot(ds, [inls[s]["obj_cent"] for s in ds], linestyle="--",marker="o", label=strdist)
        else
            plot([s for s in ds if inls[s]["obj_dec"]>0.1], [inls[s]["obj_dec"] for s in ds if inls[s]["obj_dec"]>0.1], linestyle="--",marker="o", label=strdist)
        end

    end
    legend()
    if version == "centralized"
        outstr = "plots/obj_cent_HM_"
    else
        outstr = "plots/obj_dec_HM_"
    end
    if expostir == true
        outstr = string(outstr, "expostir_", els, ".pdf")
    else
        outstr = string(outstr, els, ".pdf")
    end


    xlabel(L"$\delta$")
    ylabel("Objective")
    savefig(outstr)
    show()
end

function PlotAvgDemand(distr, types, elastic, version)
    els = nothing
    if elastic
        els = "elastic"
    else
        els = "inelastic"
    end
    for k=1:length(distr)
        dst = distr[k]
        println("========================")
        println("Distribution: ", dst)
        println("========================")
        fm = Dict(1=>dst, 2=>dst)
        strdist = join(dst,'_')
        infile = string("outputs/simulations_outcome_HM_", join(dst,'_'))
        if elastic
            els = "elastic"
            infile = string(infile, "_elastic.txt")
        else
            els = "inelastic"
            infile = string(infile, "_inelastic.txt")
        end
        sol = ReadOutputHM(infile)

        nsupp, ntypes, nvars, sts, G_x, G_t, h, A, b, q_t, f, Theta = InputsConstraintsCentralized(types, fm)

        skeys =  sort([i for i in keys(sol)])

        d1_cent = [sum([sol[skeys[r]]["x_cent"][nsupp*(s-1)+1]*f[Theta[s]] for s=1:sts ]) for r=1:length(skeys) ]
        d1_dec = [sum([sol[skeys[r]]["x_dec"][nsupp*(s-1)+1]*f[Theta[s]] for s=1:sts ]) for r=1:length(skeys) ]
        if version == "centralized"
            plot(skeys, d1_cent, linestyle="--",marker="o", label=strdist)
        else
            plot(skeys, d1_dec, linestyle="--",marker="o", label=strdist)
        end

    end
    legend()
    if version == "centralized"
        outstr = string("plots/avg_demand_cent_HM_", els, ".pdf")
    else
        outstr = string("plots/avg_demand_dec_HM_", els, ".pdf")
    end
    xlabel(L"$\delta$")
    ylabel("Objective")
    savefig(outstr)
    show()
end

function PlotVarTransfers(distr, types, elastic, version, expostir)
    els = nothing
    if elastic
        els = "elastic"
    else
        els = "inelastic"
    end
    for k=1:length(distr)
        dst = distr[k]
        println("========================")
        println("Distribution: ", dst)
        println("========================")
        fm = Dict(1=>dst, 2=>dst)
        strdist = join(dst,'_')
        infile = string("outputs/simulations_outcome_HM_", join(dst,'_'))
        if expostir
            infile = string(infile, "_expostir")
        end

        if elastic
            els = "elastic"
            infile = string(infile, "_elastic.txt")
        else
            els = "inelastic"
            infile = string(infile, "_inelastic.txt")
        end
        sol = ReadOutputHM(infile)

        nsupp, ntypes, nvars, sts, G_x, G_t, h, A, b, q_t, f, Theta = InputsConstraintsCentralized(types, fm)

        skeys =  sort([i for i in keys(sol)])
        avgT_cent = Dict(r=>mean([sol[skeys[r]]["t_cent"][nsupp*(s-1)+1] for s=1:sts ]) for r=1:length(skeys) )
        avgT_dec = Dict(r=>mean([sol[skeys[r]]["t_dec"][nsupp*(s-1)+1] for s=1:sts ]) for r=1:length(skeys) )
        d1_cent = [sum([(sol[skeys[r]]["t_cent"][nsupp*(s-1)+1]-avgT_cent[r]  )^2*f[Theta[s]] for s=1:sts ]) for r=1:length(skeys) ]
        d1_dec = [sum([(sol[skeys[r]]["t_dec"][nsupp*(s-1)+1]-avgT_dec[r]  )^2*f[Theta[s]] for s=1:sts ]) for r=1:length(skeys) ]
        if version == "centralized"
            plot(skeys, d1_cent, linestyle="--",marker="o", label=strdist)
        else
            plot(skeys, d1_dec, linestyle="--",marker="o", label=strdist)
        end

    end
    legend()
    if version == "centralized"
        outstr = "plots/variance_transfers_cent_HM_"
    else
        outstr = "plots/variance_transfers_dec_HM_"
    end
    if expostir == true
        outstr = string(outstr, "expostir_", els, ".pdf")
    else
        outstr = string(outstr, els, ".pdf")
    end
    xlabel(L"$\delta$")
    ylabel("Objective")
    savefig(outstr)
    show()
end

function PlotNumberActiveSuppliers(distr, types, elastic, version)
    els = nothing
    if elastic
        els = "elastic"
    else
        els = "inelastic"
    end
    for k=1:length(distr)
        dst = distr[k]
        println("========================")
        println("Distribution: ", dst)
        println("========================")
        fm = Dict(1=>dst, 2=>dst)
        strdist = join(dst,'_')
        infile = string("outputs/simulations_outcome_HM_", join(dst,'_'))
        if elastic
            els = "elastic"
            infile = string(infile, "_elastic.txt")
        else
            els = "inelastic"
            infile = string(infile, "_inelastic.txt")
        end
        sol = ReadOutputHM(infile)

        nsupp, ntypes, nvars, sts, G_x, G_t, h, A, b, q_t, f, Theta = InputsConstraintsCentralized(types, fm)

        skeys =  sort([i for i in keys(sol)])

        d1_cent = [sum([ sum([sol[skeys[r]]["x_cent"][nsupp*(s-1)+1]> 1e-2])*f[Theta[s]] for s=1:sts ]) for r=1:length(skeys) ]
        d1_dec = [sum([ sum([sol[skeys[r]]["x_dec"][nsupp*(s-1)+1]> 1e-2])*f[Theta[s]] for s=1:sts ]) for r=1:length(skeys) ]
        if version == "centralized"
            plot(skeys, d1_cent, linestyle="--",marker="o", label=strdist)
        else
            plot(skeys, d1_dec, linestyle="--",marker="o", label=strdist)
        end

    end
    legend()
    if version == "centralized"
        outstr = string("plots/active_suppliers_cent_HM_", els, ".pdf")
    else
        outstr = string("plots/active_suppliers_dec_HM_", els, ".pdf")
    end
    xlabel(L"$\delta$")
    ylabel("Objective")
    savefig(outstr)
    show()
end

function PlotOwnElasticities(distr, types, elastic, version)
    els = nothing
    if elastic
        els = "elastic"
    else
        els = "inelastic"
    end
    for k=1:length(distr)
        dst = distr[k]
        println("========================")
        println("Distribution: ", dst)
        println("========================")
        strdist = join(dst,'_')
        infile = string("outputs/simulations_outcome_HM_", join(dst,'_'))
        if elastic
            els = "elastic"
            infile = string(infile, "_elastic.txt")
        else
            els = "inelastic"
            infile = string(infile, "_inelastic.txt")
        end
        own_e, cross_e, own_p, cross_p = Elasticities(infile, types)
        ds = sort([i for i in keys(own_e)])
        els_cent = [own_e[s][1]["cent"]/own_p[s][1]["cent"] for s in ds]
        els_dec = [own_e[s][1]["dec"]/own_p[s][1]["dec"] for s in ds]
        println("===>", [own_e[s][1]["cent"]/own_p[s][1]["cent"] for s in ds])
        println("===>", [own_e[s][1]["dec"]/own_p[s][1]["dec"] for s in ds])
        if version == "centralized"
            plot([i for i=1:length(els_cent) if els_cent[i] > -100], [els_cent[i] for i=1:length(els_cent) if els_cent[i] > -100], linestyle="--",marker="o", label=strdist)
        else
            plot([i for i=1:length(els_dec) if els_dec[i] > -100], [els_dec[i] for i=1:length(els_dec) if els_dec[i] > -100], linestyle="--",marker="o", label=strdist)
        end

    end
    legend()
    if version == "centralized"
        outstr = string("plots/own_elast_cent_HM_", els, ".pdf")
    else
        outstr = string("plots/own_elast_dec_HM_", els, ".pdf")
    end
    xlabel(L"$\delta$")
    ylabel("Objective")
    savefig(outstr)
    show()
end

#### MAIN ####
types = Dict(1=>10, 2=>12)
distributions = [[0.1, 0.9], [0.25, 0.75], [0.4, 0.6], [0.5,0.5], [0.6,0.4], [0.75, 0.25], [0.9, 0.1]]
elastic = false
version = "decentralized"
expostir = true

# sol = ReadOutputHM("outputs/simulations_outcome_HM_0.1_0.9_inelastic.txt")
# ut, outval, outprob = CheckExPostIR("outputs/simulations_outcome_HM_0.1_0.9_inelastic.txt", types)
#
# for k in keys(ut)
#     println("================")
#     println("Delta = ", k)
#     println("================")
#     println([ut[k][s][1]["cent"] for s=1:4])
#     println([sol[k]["t_cent"][2*(s-1)+1] for s=1:4])
#     println([sol[k]["x_cent"][2*(s-1)+1] for s=1:4])
# end

PlotObjectives(distributions, elastic, version, expostir)
# PlotOwnElasticities(distributions, types, elastic, version)
# PlotAvgDemand(distributions, types, elastic, version)
# PlotNumberActiveSuppliers(distributions, types, elastic, version)
# own_e, cross_e, own_p, cross_p = Elasticities("outputs/simulations_outcome_HM_0.1_0.9_inelastic.txt", types)
#
# sol = ReadOutputHM("outputs/simulations_outcome_HM_0.1_0.9_inelastic.txt")
# for k in keys(sol)
#     println("============================")
#     println("delta = ", k)
#     println("============================")
#     println(k, " ", sol[k]["x_cent"])
#     println(k, " ", sol[k]["x_dec"])
# end



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

# own_e, cross_e, own_p, cross_p = Elasticities("outputs/simulations_outcome_HM_0.9_0.1_inelastic.txt", types)
# println(own_e)
# println(own_p)
# for k in keys(own_e)
#     for i in keys(own_e[k])
#         println(k, " ", i, " ", "cent", " ", own_e[k][i]["cent"])
#     end
# end
