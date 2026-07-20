# Convergence worker for the ORIGINAL (unoptimized) Julia RTIGER core, at an
# arbitrary rigidity. Fits one core to convergence on the shared-3 Arabidopsis panel
# decimated by odd index `level`, capturing the per-iteration progress log and a
# params/Viterbi dump (same format as the preserved r=2 runs in data/bench_ref/).
# Used to regenerate the original-Julia baseline at the r=250 operating point for the
# small sizes we then project from. Adapted from the fork's 33_panel_fit.jl.
#   julia scripts/rtiger_julia_conv_worker.jl <srcdir> <level 0..4> <rigidity>
using LinearAlgebra, Printf, DelimitedFiles
srcdir = ARGS[1]
level = parse(Int, ARGS[2])
rig = parse(Int, ARGS[3])
include(joinpath(srcdir, "AuxilaryFunctions.jl"))
include(joinpath(srcdir, "rHMM_methods.jl"))
core = occursin("orig", srcdir) ? "orig" : "opt"

# PRE-SPLIT: this worker loads ONLY the thinned panel for the level it benchmarks
# (thin_L<level>/, written by scripts/materialize_thinned_panel.R); no in-script thinning.
paneldir = "/Users/fvrodriguez/repos/zealhmm/data/rtiger_shared3_input/thin_L$(level)"
outdir = "/Users/fvrodriguez/repos/zealhmm/results/bench/orig_conv_r$(rig)"
mkpath(outdir)
files = ["S1" => "sampleBN.txt", "S2" => "sampleZ.txt", "S3" => "sampleAU.txt"]
chrs = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
eps = 0.01
max_iter = 50

raw = Dict{String,Matrix{Any}}()
for (s, f) in files
    raw[s] = readdlm(joinpath(paneldir, f))
end
N = size(raw["S1"], 1)
chrcol = string.(raw["S1"][:, 1])

function rinit(r)
    a = fill(0.1, 3, 3) + 10 * Matrix{Float64}(I, 3, 3)
    a = a ./ sum(a, dims = 2)
    Dict(:logpi => log.([1 / 3, 1 / 3, 1 / 3]), :transition => a, :logtransition => log.(a),
        :paraBetaAlpha => [20.0, 20.0, 1.0], :paraBetaBeta => [1.0, 20.0, 20.0],
        :nstates => 3, :rigidity => r)
end
function build(ix)
    O = Dict{Any,Any}()
    for (s, _) in files
        M = raw[s]
        O[s] = Dict{Any,Any}()
        for c in chrs
            rows = ix[chrcol[ix].==c]
            isempty(rows) && continue
            k = Int.(round.(M[rows, 4]))
            n = k .+ Int.(round.(M[rows, 6]))
            O[s][c] = hcat(k, n)
        end
    end
    O
end

fit(build(collect(1:600)), nothing, rinit(rig), 2, eps, false, true, false, 20, nothing, true)  # warmup (>=2r markers)

O = build(collect(1:N))   # the FULL pre-split panel for this level (no thinning)
mps = N
logpath = joinpath(outdir, "log_$(core)_conv_L$(level).log")
isfile(logpath) && rm(logpath)
GC.gc()
t0 = time()
res = fit(O, nothing, rinit(rig), max_iter, eps, false, true, false, 20, nothing, true;
    progress_log = logpath)
rt = time() - t0
it = res[:numberofiterations]

open(joinpath(outdir, "conv_$(core)_L$(level).txt"), "w") do f
    p = res[:parameterSet]
    println(f, "level=", level, " mps=", mps, " iters=", it, " runtime=", round(rt, digits = 3))
    println(f, "alpha=", round.(vec(p[:paraBetaAlpha]); digits = 6))
    println(f, "beta=", round.(vec(p[:paraBetaBeta]); digits = 6))
    println(f, "pi=", round.(vec(p[:pi]); digits = 6))
    println(f, "transition=", round.(vec(p[:transition]); digits = 6))
    for s in sort(collect(keys(O))), c in sort(collect(keys(O[s])))
        println(f, "vit[$s][$c]=", res[:viterbiPath][s][c])
    end
end
@printf("SWEEP core=%s level=%d mps=%d iters=%d runtime=%.3f rig=%d\n", core, level, mps, it, rt, rig)
