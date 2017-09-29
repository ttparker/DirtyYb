const dir = "JNeg/Delta0p4/L46/"  # directory storing the raw data files
const realizations = ["Real1", "Real2", "Real3", "Real4", "Real5"]  # Raw data output files for each disorder realization

import JLD
include("CustomTypes.jl")

mkpath("Analyzed/" * dir)
for realization in realizations
  system = SystemSummary(JLD.load("Results/" * dir * realization * ".jld", "params")::SystemParameters)
  for run in JLD.load("Results/" * dir * realization * ".jld", "runs")::Vector{RawRunData}
    push!(system.Ts, run.T)
    push!(system.avgEnergies, mean(run.energyList))
    push!(system.energyVariances, var(run.energyList))
    push!(system.absAvgPsis, abs(mean(run.psiList)))
    push!(system.avgAbsPsis, mean(abs(run.psiList)))
    push!(system.absAvgPsi2s, abs(mean(run.psi2List)))
    push!(system.avgAbsPsi2s, mean(abs(run.psi2List)))
    push!(system.absAvgPsi3s, abs(mean(run.psi3List)))
    push!(system.avgAbsPsi3s, mean(abs(run.psi3List)))
  end
  JLD.save("Analyzed/" * dir * realization * ".jld", "system", system)
end
