const dataFiles = ["Output"]  # Raw data output files for this system to be read
const filename = "L46"  # Name of file to save this system's data
const legendLabel = "N = 2116"  # How to label this system's data in plots

import JLD
include("CustomTypes.jl")

function analyzeRun!(f::IOStream, system::SystemSummary)
  line = readline(f)
  if line == "Summary:\n"
    return true
  else
    T = parse(Float64, split(line)[2])
    push!(system.Ts, T)
    readline(f)  # "Energies:"
    energies = map(x -> parse(Float64, x), split(readline(f), ", "))
    push!(system.avgEnergies, mean(energies))
    push!(system.energyVariances, var(energies))
    readline(f)  # "Psis:"
    psis = eval.(parse.(split(readline(f), ", ")))
    push!(system.absAvgPsis, abs(mean(psis)))
    push!(system.avgAbsPsis, mean(abs(psis)))
    push!(system.absAvgPsi2s, abs(mean(psis.^2)))
    push!(system.avgAbsPsi2s, mean(abs2(psis)))
    readline(f)
    println("Run T = ", T, " analyzed")
    return false
  end
end

system = SystemSummary(legendLabel)
for dataFile in dataFiles
  f = open(dataFile, "r")
  summarySection = false
  while !eof(f) && !summarySection
    summarySection = analyzeRun!(f, system)
  end
  if dataFile == dataFiles[end] && summarySection  # use summary section from last file to fill out system.params
    hamParams = split(readline(f), r": |, ")[[2, 4, 6, 8, 10, 12]]
    system.params.Jzz, system.params.Jpm, system.params.Jpmpm, system.params.Jzpm , system.params.Delta = map(x -> parse(Float64, x), hamParams[1:5])
    system.params.L = parse(Int64, hamParams[6])
    system.params.thermalizationSweeps, system.params.equilibriumSweeps = map(x -> parse(Int64, x), split(readline(f), r": |, ")[[2, 4]])
  end
  close(f)
end
mkpath("Analyzed")
JLD.save("Analyzed/" * filename * ".jld", filename, system)
