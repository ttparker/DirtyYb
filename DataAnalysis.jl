const dir = "Delta0p4/L46/"  # directory storing the raw data files
const realizations = ["Real1", "Real2", "Real3", "Real4", "Real5", "Real6"]  # Raw data output files for each disorder realization
const legendLabel = "N = 2116"  # How to label this system's data in plots

import JLD
include("CustomTypes.jl")

for realization in realizations
  system = SystemSummary(legendLabel)
  system.params = JLD.load(dir * realization * ".jld", "params")::SystemParameters
  for run in JLD.load(dir * realization * ".jld", "runs")::Vector{RawRunData}
    push!(system.Ts, run.T)
    push!(system.avgEnergies, mean(run.energyList))
    push!(system.energyVariances, var(run.energyList))
    push!(system.absAvgPsis, abs(mean(run.psiList)))
    push!(system.avgAbsPsis, mean(abs(run.psiList)))
    push!(system.absAvgPsi2s, abs(mean(run.psiList.^2)))
    push!(system.avgAbsPsi2s, mean(abs2(run.psiList)))
  end
  mkpath("Analyzed/" * dir * realization * "/")
  JLD.save("Analyzed/" * dir * realization * "/" * realization * ".jld", "system", system)
end
