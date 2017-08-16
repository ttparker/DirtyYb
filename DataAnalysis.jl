const dataFiles = ["Output.jld"]  # Raw data output files for this system to be read
const filename = "L46"  # Name of file to save this system's data
const legendLabel = "N = 2116"  # How to label this system's data in plots

import JLD
include("CustomTypes.jl")

system = SystemSummary(legendLabel)
for dataFile in dataFiles
  runs = JLD.load(dataFile)["runs"]
  for run in runs
    push!(system.Ts, run.T)
    push!(system.avgEnergies, mean(run.energyList))
    push!(system.energyVariances, var(run.energyList))
    push!(system.absAvgPsis, abs(mean(run.psiList)))
    push!(system.avgAbsPsis, mean(abs(run.psiList)))
    push!(system.absAvgPsi2s, abs(mean(run.psiList.^2)))
    push!(system.avgAbsPsi2s, mean(abs2(run.psiList)))
  end
end
system.params = JLD.load(dataFiles[1])["params"]

mkpath("Analyzed")
JLD.save("Analyzed/" * filename * ".jld", filename, system)
