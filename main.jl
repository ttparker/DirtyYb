include("CustomTypes.jl")  # non-built-in Julia types
include("Yb.jl")  # details of this specific system
include("MC.jl")  # MC algorithm

# Read in run parameters from "Input":
params = SystemParameters()
f = open("Input", "r")
params.Jzz, params.Jpm, params.Jpmpm, params.Jzpm = [parse(Float64, s) for s in split(readline(f))]
params.L = parse(Int64, readline(f))
Ts = [parse(Float64, s) for s in split(readline(f))]
readline(f)
params.thermalizationSweeps, params.equilibriumSweeps = [eval(parse(s)) for s in split(readline(f))]
close(f)

# Read in disorder realization:
import JLD
realization = JLD.load("Realization.jld")
if realization["L"] == params.L
  params.Delta = realization["Delta"]
  bondMultipliers = realization["multipliers"]
else
  println("System size mismatch")
  quit()
end

f = open("Output","w")
for currentT in Ts  # MC runs
  params.T = currentT
  MCRun(params, bondMultipliers, f)
  println("MC run complete")
end
println(f, "Summary:\nJzz: ", params.Jzz, ", Jpm: ", params.Jpm, ", Jpmpm: ", params.Jpmpm, ", Jzpm: ", params.Jzpm, ", Delta: ", params.Delta, ", L: ", params.L, "\nThermalization sweeps: ", params.thermalizationSweeps, ", Equilibrium sweeps: ", params.equilibriumSweeps)
close(f)
println("Done")
