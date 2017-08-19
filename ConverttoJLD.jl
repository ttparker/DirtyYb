const dataFiles = ["Output"]  # Raw data output files for this system to be read
const filename = "L10"  # Name of file to save this system's data
const legendLabel = "N = 100"  # How to label this system's data in plots

import JLD
include("CustomTypes.jl")

function analyzeRun!(f::IOStream, runs::Vector{RawRunData})
  line = readline(f)
  if line == "Summary:\n"
    return true
  else
    T = parse(Float64, split(line)[2])
    readline(f)  # "Energies:"
    energies = map(x -> parse(Float64, x), split(readline(f), ", "))
    readline(f)  # "Psis:"
    psis = eval.(parse.(split(readline(f), ", ")))
    readline(f)
    push!(runs, RawRunData(T, energies, psis))
    println("Run T = ", T, " read")
    return false
  end
end

for dataFile in dataFiles
  f = open(dataFile, "r")
  summarySection = false
  runs = Vector{RawRunData}()
  while !eof(f) && !summarySection
    summarySection = analyzeRun!(f, runs)
  end
  params = SystemParameters()
  hamParams = split(readline(f), r": |, |\n")[[2, 4, 6, 8, 10, 12]]
  params.Jzz, params.Jpm, params.Jpmpm, params.Jzpm, params.Delta = map(x -> parse(Float64, x), hamParams[1:5])
  params.L = parse(Int64, hamParams[6])
  params.thermalizationSweeps, params.equilibriumSweeps = map(x -> parse(Int64, x), split(readline(f), r": |, |\n")[[2, 4]])
  close(f)
  JLD.save(dataFile * ".jld", "params", params, "runs", runs)
  println(dataFile, " analyzed")
end
println("Done")
