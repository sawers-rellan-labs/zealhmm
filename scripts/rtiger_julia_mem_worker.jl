# Peak-RSS worker: original (array-retaining) Julia RTIGER core, one EM iteration
# (max_iter=1 -> the E-step allocates all per-sample arrays, so peak RSS is reached)
# on the Arabidopsis shared-3 panel thinned by odd index `level`. Peak RSS is
# measured by the parent (/usr/bin/time -l). Usage: julia rtiger_julia_mem_worker.jl L
using DelimitedFiles, LinearAlgebra
SRC   = "/Users/fvrodriguez/repos/rtiger-fork-assets/agent/scale_check/orig_base"
PANEL = "/Users/fvrodriguez/repos/zealhmm/data/rtiger_shared3_input"
include(joinpath(SRC, "AuxilaryFunctions.jl"))
include(joinpath(SRC, "rHMM_methods.jl"))

level = parse(Int, ARGS[1])
files = ["S1" => "sampleBN.txt", "S2" => "sampleZ.txt", "S3" => "sampleAU.txt"]
chrs  = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
eps = 0.01; rig = 2

raw = Dict{String,Matrix{Any}}()
for (s, f) in files; raw[s] = readdlm(joinpath(PANEL, f)); end
N = size(raw["S1"], 1); chrcol = string.(raw["S1"][:, 1])
idx = let ix = collect(1:N)
    for _ in 1:level; ix = ix[1:2:end]; end
    ix
end

function rinit(r)
    a = fill(0.1, 3, 3) + 10 * Matrix{Float64}(I, 3, 3); a = a ./ sum(a, dims = 2)
    Dict(:logpi => log.([1/3, 1/3, 1/3]), :transition => a, :logtransition => log.(a),
         :paraBetaAlpha => [20.0, 20.0, 1.0], :paraBetaBeta => [1.0, 20.0, 20.0],
         :nstates => 3, :rigidity => r)
end
function build(ix)
    O = Dict{Any,Any}()
    for (s, _) in files
        M = raw[s]; O[s] = Dict{Any,Any}()
        for c in chrs
            rows = ix[chrcol[ix] .== c]; isempty(rows) && continue
            k = Int.(round.(M[rows, 4])); n = k .+ Int.(round.(M[rows, 6]))
            O[s][c] = hcat(k, n)
        end
    end
    O
end

O = build(idx); mps = length(idx)
# fit(Obs, info, init, max_iter, eps, trace, all, random, nsamples, specific, post_processing, DEBUG)
fit(O, nothing, rinit(rig), 1, eps, false, true, false, 20, nothing, true, false)
println("MEMJUL_DONE level=$level markers=$mps")
