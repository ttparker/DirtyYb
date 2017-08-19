const reals = ["Real1", "Real2", "Real3", "Real4", "Real5", "Real6"]
const blocks = ["Block1", "Block2", "Block3"]

import JLD
include(reals[1] * "/CustomTypes.jl")

for real in reals
  dic = JLD.load(real * "/" * blocks[1] * ".jld")
  params = dic["params"]::SystemParameters
  combinedRuns = dic["runs"]::Vector{RawRunData}
  println(blocks[1], " added")
  for block in blocks[2:end]
    append!(combinedRuns, JLD.load(real * "/" * block * ".jld")["runs"])
    println(block, " added")
  end
  JLD.save(real * ".jld", "params", params, "runs", combinedRuns)
  println(real " converted")
end
