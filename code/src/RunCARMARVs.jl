using CARMA
using CSV
using Ensemble
using HDF5
using Printf
using PyCall
using PyPlot
using Statistics
using TOML

@pyimport seaborn as sns
@pyimport corner

sns.set_context("notebook")
sns.set_style("ticks")
sns.set_palette("colorblind")

function parse_optfile(option_file)
    opts = TOML.parsefile(option_file)
    rv_opts = opts["RV"]
    p_opts = opts["Planet"]
    c_opts = opts["CARMA"]
    o_opts = opts["Output"]
    s_opts = opts["Sampler"]

    (rv_opts, p_opts, c_opts, o_opts, s_opts)
end

function read_rvs(rv_dir, tkey, rvkey, drvkey; delim=',', ignorerepeated=false)
    ts = Array{Float64, 1}[]
    vs = Array{Float64, 1}[]
    dvs = Array{Float64, 1}[]
    for f in readdir(rv_dir)
        d = CSV.read(joinpath(rv_dir, f), delim=delim, ignorerepeated=ignorerepeated)

        push!(ts, d[tkey])
        push!(vs, d[rvkey])
        push!(dvs, d[drvkey])
    end

    return (ts, vs, dvs)
end

function make_posterior(rv_opts, p_opts, c_opts)
    rv_dir = rv_opts["rvdir"]

    delim = get(rv_opts, "delim", ',')
    ignorerepeated = get(rv_opts, "ignorerepeated", false)

    ts, vs, dvs = read_rvs(rv_dir, Symbol(rv_opts["timekey"]), Symbol(rv_opts["rvkey"]), Symbol(rv_opts["rverrkey"]), delim=delim, ignorerepeated=ignorerepeated)

    post = CARMAKepler.MultiEpochPosterior(ts, vs, dvs, p_opts["Pmin"], p_opts["Pmax"], c_opts["ndrw"], c_opts["nosc"], c_opts["rmin"], c_opts["rmax"], c_opts["fmin"], c_opts["fmax"], c_opts["Qmax"])

    post
end

function load_run(option_file)
    rv_opts, p_opts, c_opts, o_opts, s_opts = parse_optfile(option_file)
    post = make_posterior(rv_opts, p_opts, c_opts)

    h5open(joinpath(o_opts["outdir"], o_opts["statefile"]), "r") do f
        global ps, ns = read(f, post)
    end

    (post, ps, ns)
end

function run_carma_rvs(option_file)
    rv_opts, p_opts, c_opts, o_opts, s_opts = parse_optfile(option_file)

    mkpath(o_opts["outdir"])
    cfile = joinpath(o_opts["outdir"], o_opts["ckptfile"])
    ofile = joinpath(o_opts["outdir"], o_opts["statefile"])

    post = make_posterior(rv_opts, p_opts, c_opts)

    if isfile(cfile)
        h5open(cfile, "r") do f
            ns = EnsembleNest.NestState(f, logl=x->CARMAKepler.log_likelihood(post, x), logp=x->CARMAKepler.log_prior(post, x))
        end
    else
        ns = EnsembleNest.NestState(x->CARMAKepler.log_likelihood(post, x), x->CARMAKepler.log_prior(post, x), CARMAKepler.draw_prior(post, s_opts["nlive"]), s_opts["nmcmc"])
    end

    EnsembleNest.run!(ns, s_opts["dZstop"], verbose=true, ckpt_file=cfile)

    h5open(joinpath(o_opts["outdir"], o_opts["statefile"]), "w") do f
        write(f, post, ns)
    end
    rm(cfile)

    post
end

function make_default_plots(option_file)
    rv_opts, p_opts, c_opts, o_opts, s_opts = parse_optfile(option_file)

    sfile = joinpath(o_opts["outdir"], o_opts["statefile"])

    post, ps, ns = load_run(option_file)

    # CARMA Parameters
    drw_rms = [p.drw_rms for p in ps]
    if size(drw_rms[1], 1) == 1
        drw_rms = vcat(drw_rms...)
        drw_rms_label = [L"$A_\mathrm{drw}$"]
    else
        drw_rms_label = [@sprintf("\$A_\\mathrm{drw}^{(%d)}\$", i) for i in 1:size(drw_rms[1],1)]
    end

    r_drw = [p.drw_rate for p in ps]
    if size(r_drw[1], 1) == 1
        r_drw = vcat(r_drw...)
        r_drw_label = [L"$r_\mathrm{drw}$"]
    else
        r_drw_label = [@sprintf("\$r_\\mathrm{drw}^{(%d)}\$", i) for i in 1:size(r_drw[1], 1)]
    end

    osc_rms = [p.osc_rms for p in ps]
    if size(osc_rms[1],1) == 1
        osc_rms = vcat(osc_rms...)
        osc_rms_label = [L"$A_\mathrm{osc}$"]
    else
        osc_rms_label = [@sprintf("\$A_\\mathrm{osc}^{(%d)}\$", i) for i in 1:size(osc_rms[1],1)]
    end

    osc_freq = [p.osc_freq for p in ps]
    if size(osc_freq[1],1) == 1
        osc_freq = vcat(osc_freq...)
        osc_freq_label = [L"$f$"]
    else
        osc_freq_label = [@sprintf("\$f^{(%d)}\$", i) for i in 1:size(osc_freq[1],1)]
    end

    osc_Q = [p.osc_Q for p in ps]
    if size(osc_Q[1],1) == 1
        osc_Q = vcat(osc_Q...)
        osc_Q_label = [L"$Q$"]
    else
        osc_Q_label = [@sprintf("\$Q^{(%d)}\$", i) for i in 1:size(osc_Q[1],1)]
    end

    d = hcat(drw_rms, r_drw, osc_rms, osc_freq, osc_Q)
    l = vcat(drw_rms_label, r_drw_label, osc_rms_label, osc_freq_label, osc_Q_label)
    corner.corner(d, labels=l)

    tight_layout()
    savefig(joinpath(o_opts["outdir"], "CARMA-params.pdf"))

    # Kepler Params
    corner.corner(hcat([p.K for p in ps], [p.P for p in ps], [p.e for p in ps]),
                  labels=[L"$K$", L"$P$", L"$e$"])
    tight_layout()
    savefig(joinpath(o_opts["outdir"], "Kepler-params.pdf"))

    # nu Figure
    figure()
    d = hcat([p.nu for p in ps]...)
    for j in 1:size(d, 1)
        sns.distplot(vec(d[j, :]))
    end
    xlabel(L"$\nu$")
    ylabel(L"$p\left( \nu \right)$")
    tight_layout()
    savefig(joinpath(o_opts["outdir"], "nu.pdf"))

    # Three draws from the posterior over the RV model
    tmin = minimum([minimum(t) for t in post.ts])
    tmax = maximum([maximum(t) for t in post.ts])
    tpred = collect(range(tmin, stop=tmax, length=1000))
    for draw in 1:3
        i = rand(1:size(ps,1))
        vp, vvp = CARMAKepler.predict(post, ps[i], tpred)

        figure()
        for (j, (t, v, dv)) in enumerate(zip(post.ts, post.ys, post.dys))
            errorbar(t, v .- ps[i].mu[j], dv, fmt=".")
        end

        plot(tpred, vp, "-k")
        fill_between(tpred, vp.+sqrt.(vvp), vp.-sqrt.(vvp), color="k", alpha=0.25)
        fill_between(tpred, vp.+ 2 .* sqrt.(vvp), vp.- 2 .* sqrt.(vvp), color="k", alpha=0.25)

        xlabel(L"$t$")
        ylabel(L"$v_r$")
        tight_layout()
        savefig(joinpath(o_opts["outdir"], @sprintf("draw-%d.pdf", draw)))
    end

    # And the PSD
    T = tmax - tmin
    df = 1.0/T
    dt_min = minimum(diff(sort(vcat(post.ts...))))
    fmax = 1.0/(2.0*dt_min)

    fs = exp.(range(log(df/10.0), stop=log(fmax*10.0), length=1000))
    psds = [CARMAKepler.psd(post, p, fs) for p in ps]

    m = zeros(1000)
    l = zeros(1000)
    h = zeros(1000)
    ll = zeros(1000)
    hh = zeros(1000)

    for i in 1:1000
        v = [p[i] for p in psds]
        m[i] = median(v)
        l[i] = quantile(v, 0.16)
        h[i] = quantile(v, 0.84)
        ll[i] = quantile(v, 0.025)
        hh[i] = quantile(v, 0.975)
    end

    figure()
    loglog(fs, m)
    fill_between(fs, h, l, color=sns.color_palette()[1], alpha=0.25)
    fill_between(fs, hh, ll, color=sns.color_palette()[1], alpha=0.25)

    xlabel(L"$f$")
    ylabel(L"$P(f)$")

    tight_layout()
    savefig(joinpath(o_opts["outdir"], "psd.pdf"))
end
