const dirs = ["Light/JNeg/Delta0p3/L46/", "Light/JNeg/Delta0p3/L64/", "Light/JNeg/Delta0p3/L90/", "Light/JNeg/Delta0p3/L128/"#=, "Light/JPos/Delta0p2/L46/", "Light/JPos/Delta0p2/L64/", "Light/JPos/Delta0p2/L90/", "Light/JPos/Delta0p2/L128/"=#]
const realizationNames = ["Real1", "Real2", "Real3", "Real4", "Real5", "Real6", "Real7", "Real8", "Real9", "Real10"]
const ConvergenceThreshold = .2

import JLD
include("CustomTypes.jl")

for dir in dirs
  # load the realizations:
  real1 = JLD.load("Analyzed/" * dir * realizationNames[1] * ".jld", "system")::SystemSummary
  Ts = real1.Ts
  reals = [real1]
  for i in 2:length(realizationNames)
    real = JLD.load("Analyzed/" * dir * realizationNames[i] * ".jld", "system")::SystemSummary
    if real.Ts == Ts
      push!(reals, real)
    else
      println("Error: real ", i, "'s Ts != first run's Ts")
      quit()
    end
  end

  # find the converged runs for each temperature:
  averaged = SystemSummary(real1.params)
  averaged.Ts = Ts
  for TNo in 1:length(Ts)
    Cs = [real.energyVariances[TNo] * averaged.params.L^2 / Ts[TNo]^2 for real in reals]
    convergedRealNs = find(C -> 1 - ConvergenceThreshold <= C / mean(Cs) <= 1 + ConvergenceThreshold, Cs)
    print(length(convergedRealNs), " ")
    push!(averaged.avgEnergies, mean([real.avgEnergies[TNo] for real in reals[convergedRealNs]]))
    push!(averaged.energyVariances, mean([real.energyVariances[TNo] for real in reals[convergedRealNs]]))
    push!(averaged.absAvgPsis, mean([real.absAvgPsis[TNo] for real in reals[convergedRealNs]]))
    push!(averaged.avgAbsPsis, mean([real.avgAbsPsis[TNo] for real in reals[convergedRealNs]]))
    push!(averaged.absAvgPsi2s, mean([real.absAvgPsi2s[TNo] for real in reals[convergedRealNs]]))
    push!(averaged.avgAbsPsi2s, mean([real.avgAbsPsi2s[TNo] for real in reals[convergedRealNs]]))
    push!(averaged.absAvgPsi3s, mean([real.absAvgPsi3s[TNo] for real in reals[convergedRealNs]]))
    push!(averaged.avgAbsPsi3s, mean([real.avgAbsPsi3s[TNo] for real in reals[convergedRealNs]]))
  end
  println()
  JLD.save("Analyzed/" * dir * "ConvergedAveraged.jld", "system", averaged)
end
