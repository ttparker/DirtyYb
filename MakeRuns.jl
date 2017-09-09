const L = 46
const Delta = 0.4
const Jpmpm = -0.9
const runLabels = ["Real1", "Real2", "Real3", "Real4", "Real5"]

if L == 46
  const blocks = [1.00:.01:1.24, 1.25:.01:1.49, 1.50:.01:1.74, 1.75:.01:2.00]
elseif L == 64
  const blocks = [1.00:.01:1.12, 1.13:.01:1.24, 1.25:.01:1.36, 1.37:.01:1.49, 1.50:.01:1.62, 1.63:.01:1.74, 1.75:.01:1.87, 1.88:.01:2.00]
elseif L == 90
  const blocks = [1.00:.01:1.06, 1.07:.01:1.12, 1.13:.01:1.18, 1.19:.01:1.24, 1.25:.01:1.30, 1.31:.01:1.36, 1.37:.01:1.43, 1.44:.01:1.49, 1.50:.01:1.56, 1.57:.01:1.62, 1.63:.01:1.68, 1.69:.01:1.74, 1.75:.01:1.81, 1.82:.01:1.87, 1.88:.01:1.94, 1.95:.01:2.00]
elseif L == 128
  const blocks = [1.00:.01:1.02, 1.03:.01:1.05, 1.06:.01:1.08, 1.09:.01:1.11, 1.12:.01:1.14, 1.15:.01:1.17, 1.18:.01:1.20, 1.21:.01:1.23, 1.24:.01:1.26, 1.27:.01:1.29, 1.30:.01:1.32, 1.33:.01:1.35, 1.36:.01:1.38, 1.39:.01:1.41, 1.42:.01:1.44, 1.45:.01:1.47, 1.48:.01:1.50, 1.51:.01:1.53, 1.54:.01:1.56, 1.57:.01:1.59, 1.60:.01:1.62, 1.63:.01:1.65, 1.66:.01:1.68, 1.69:.01:1.71, 1.72:.01:1.74, 1.75:.01:1.77, 1.78:.01:1.80, 1.81:.01:1.83, 1.84:.01:1.86, 1.87:.01:1.89, 1.90:.01:1.92, 1.93:.01:1.95, 1.96:.01:1.98, 1.99:.01:2.0]
end

const Jzz = 1.
const Jpm = 0.915
const Jzpm = 0.5
const thermalizationSweeps = 4 * 10^5
const equilibriumSweeps = 1 * 10^6

const baseDirectory = pwd() * "/"

f = open("MakeRuns.sh", "w")
println(f, "#!/bin/bash")
println(f)
for runLabel in runLabels
  runDirectory = "Runs/" * (Jpmpm < 0 ? "JNeg/" : "JPos/") * "Delta" * replace(string(Delta), ".", "p") * "/L" * string(L) * "/" * runLabel * "/"
  if isdir(runDirectory)
    println("Skipping ", runLabel, ", name already used")
    continue
  else
    mkpath(runDirectory)
  end
  println(f, "julia GenerateDisorder.jl ", L, " ", Delta, "\n")
  for blockNo in 1:length(blocks)
    blockDirectory = runDirectory * "Block" * string(blockNo) * "/"
    println(f, "mkdir ", blockDirectory)
    println(f, "cp Realization.jld CustomTypes.jl main.jl MC.jl Yb.jl ", blockDirectory)
    println(f, "cd ", blockDirectory)

    # make Input file:
    print(f, "printf '", Jzz, " ", Jpm, " ", Jpmpm, " ", Jzpm, "\\n", L, "\\n")
    for T in blocks[blockNo]
      print(f, T, " ")
    end
    println(f, "\\n\\n", thermalizationSweeps, " ", equilibriumSweeps, "\\n' > Input")

    println(f, "printf '#!/bin/bash\\n#PBS -l nodes=1:ppn=1\\n#PBS -N ", (Jpmpm < 0 ? "IP" : "OP"), "D", Int64(10 * Delta), "L", L, runLabel, "B", blockNo, "\\n#PBS -m abe\\n#PBS -M tparker@physics.ucsb.edu\\n#PBS -j oe\\n\\ncd \$PBS_O_WORKDIR\\njulia main.jl' > job")
    println(f, "qsub job")
    println(f, "cd ", baseDirectory, "\n")
  end
  println(f, "rm Realization.jld\n")
end
close(f)
