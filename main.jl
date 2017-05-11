include("Yb.jl")  # details of this specific system
include("MC.jl")  # MC algorithm

f = open("Input", "r")
lines = readlines(f)  # read in run parameters
close(f)
currentJzz, currentJpm, currentJpmpm, currentJzpm, Delta = [parse(Float64, s) for s in split(lines[1])]
currentL = parse(Int64, lines[2])
Ts = [parse(Float64, s) for s in split(lines[3])]
thermalizationSweeps, equilibriumSweeps = [eval(parse(s)) for s in split(lines[5])]

bondMultipliers = ones(currentL, currentL, 3) + Delta * randn(currentL, currentL, 3)

f = open("Output","w")
println(f)
for currentT in Ts  # MC runs
  MCRun(currentJzz, currentJpm, currentJpmpm, currentJzpm, bondMultipliers, currentT, currentL, thermalizationSweeps, equilibriumSweeps, f)
  println("MC run complete")
end
println(f, "Jzz: ", currentJzz, ", Jpm: ", currentJpm, ", Jpmpm: ", currentJpmpm, ", Jzpm: ", currentJzpm, ", Delta: ", Delta, ", L: ", currentL, "\nBond multipliers:\n", reshape(bondMultipliers, currentL^2, 3), "\nThermalization sweeps: ", thermalizationSweeps, ", Equilibrium sweeps: ", equilibriumSweeps)
close(f)
println("Done")
