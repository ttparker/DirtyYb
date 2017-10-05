const nBlocks = 4  # L46: 4, L64: 8, L90: 16, L128: 34
const dir = "JNeg/Delta0p4/L46/"
const reals = ["Real1", "Real2", "Real3", "Real4", "Real5"]

include("CustomTypes.jl")
import JLD

mkpath("Results/" * dir)
for real in reals
  params = JLD.load("Runs/" * dir * real * "/Block1/Output.jld", "params")::SystemParameters
  runs = JLD.load("Runs/" * dir * real * "/Block1/Output.jld", "runs")::Vector{RawRunData}
  for blockNo = 2:nBlocks
    append!(runs, JLD.load("Runs/" * dir * real * "/Block" * string(blockNo) * "/Output.jld", "runs")::Vector{RawRunData})
  end
  JLD.save("Results/" * dir * real * ".jld", "params", params, "runs", runs)
end
