const dir = "Delta0p4/L46/"
const realizations = ["Real1", "Real2", "Real3", "Real4", "Real5", "Real6"]
const fileName = "L46"  # name of saved disorder-averaged data

import JLD
include("CustomTypes.jl")

disorderAveraged = JLD.load("Analyzed/" * dir * realizations[1] * "/" * realizations[1] * ".jld", "system")::SystemSummary
for realization in realizations[2:end]
  system = JLD.load("Analyzed/" * dir * realization * "/" * realization * ".jld", "system")::SystemSummary
  disorderAveraged.avgEnergies += system.avgEnergies
  disorderAveraged.energyVariances += system.energyVariances
  disorderAveraged.absAvgPsis += system.absAvgPsis
  disorderAveraged.avgAbsPsis += system.avgAbsPsis
  disorderAveraged.absAvgPsi2s += system.absAvgPsi2s
  disorderAveraged.avgAbsPsi2s += system.avgAbsPsi2s
end
const nRealizations = length(realizations)
disorderAveraged.avgEnergies /= nRealizations
disorderAveraged.energyVariances /= nRealizations
disorderAveraged.absAvgPsis /= nRealizations
disorderAveraged.avgAbsPsis /= nRealizations
disorderAveraged.absAvgPsi2s /= nRealizations
disorderAveraged.avgAbsPsi2s /= nRealizations
mkpath("Analyzed/" * dir * "Averaged/")
JLD.save("Analyzed/" * dir * "Averaged/" * fileName * ".jld", fileName, disorderAveraged)
