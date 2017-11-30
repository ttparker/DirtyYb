const dirs = ["Light/JNeg/Delta0p3/L46/", "Light/JNeg/Delta0p3/L64/", "Light/JNeg/Delta0p3/L90/", "Light/JNeg/Delta0p3/L128/", "Light/JPos/Delta0p2/L46/", "Light/JPos/Delta0p2/L64/", "Light/JPos/Delta0p2/L90/", "Light/JPos/Delta0p2/L128/"]
const realizations = ["Real1", "Real2", "Real3", "Real4", "Real5", "Real6", "Real7", "Real8", "Real9", "Real10"]

import JLD
include("CustomTypes.jl")

for dir in dirs
  disorderAveraged = JLD.load("Analyzed/" * dir * realizations[1] * ".jld", "system")::SystemSummary
  for realization in realizations[2:end]
    system = JLD.load("Analyzed/" * dir * realization * ".jld", "system")::SystemSummary
    disorderAveraged.avgEnergies += system.avgEnergies
    disorderAveraged.energyVariances += system.energyVariances
    disorderAveraged.absAvgPsis += system.absAvgPsis
    disorderAveraged.avgAbsPsis += system.avgAbsPsis
    disorderAveraged.absAvgPsi2s += system.absAvgPsi2s
    disorderAveraged.avgAbsPsi2s += system.avgAbsPsi2s
    disorderAveraged.absAvgPsi3s += system.absAvgPsi3s
    disorderAveraged.avgAbsPsi3s += system.avgAbsPsi3s
  end
  const nRealizations = length(realizations)
  disorderAveraged.avgEnergies /= nRealizations
  disorderAveraged.energyVariances /= nRealizations
  disorderAveraged.absAvgPsis /= nRealizations
  disorderAveraged.avgAbsPsis /= nRealizations
  disorderAveraged.absAvgPsi2s /= nRealizations
  disorderAveraged.avgAbsPsi2s /= nRealizations
  disorderAveraged.absAvgPsi3s /= nRealizations
  disorderAveraged.avgAbsPsi3s /= nRealizations
  JLD.save("Analyzed/" * dir * "Averaged.jld", "system", disorderAveraged, "nRealizations", nRealizations)
end
