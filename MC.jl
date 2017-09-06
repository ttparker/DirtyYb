function randomSpin(r1::Float64, r2::Float64)
  z = 2 * r1 - 1
  s = sqrt(1 - z^2)
  phi = 2 * pi * r2
  [s * cos(phi), s * sin(phi), z]
end

function MCSweep!(sys::Yb, measure::Bool)
  ns = rand(1:sys.L, sys.N, 2) # generating all rands at once saves RNG call overhead
  rands = rand(sys.N, 3)  # (for z, for phi, for flipping prob)
  for step = 1:sys.N
    n2, n3 = ns[step, :]  # pick a spin to flip
    nw = mod1(n2 + 1, sys.L)
    se = mod1(n2 - 1, sys.L)
    sw = mod1(n3 + 1, sys.L)
    ne = mod1(n3 - 1, sys.L)
    candidate = randomSpin(rands[step, 1], rands[step, 2])  # candidate new orientation
    DeltaS = candidate - sys.spins[n2, n3, :]
    DeltaE = dot(DeltaS, sys.Ja1 * (sys.bondMultipliers[se, ne, 1] * sys.spins[se, ne, :] + sys.bondMultipliers[n2, n3, 1] * sys.spins[nw, sw, :]) + sys.Ja2 * (sys.bondMultipliers[n2, n3, 2] * sys.spins[nw, n3, :] + sys.bondMultipliers[se, n3, 2] * sys.spins[se, n3, :]) + sys.Ja3 * (sys.bondMultipliers[n2, n3, 3] * sys.spins[n2, sw, :] + sys.bondMultipliers[n2, ne, 3] * sys.spins[n2, ne, :]))
    if DeltaE <= 0 || rands[step, 3] < exp(-DeltaE / sys.T)
      if measure  # update the appropriate sublattice magnetization
        n2Even = iseven(n2)
        n3Even = iseven(n3)
        if n2Even && n3Even
          sys.MA += DeltaS
        elseif !n2Even && !n3Even
          sys.MB += DeltaS
        elseif !n2Even && n3Even
          sys.MC += DeltaS
        else
          sys.MD += DeltaS
        end
      end
      sys.spins[n2, n3, :] = candidate  # replace spin by candidate
      sys.energy += DeltaE
    end
  end
end

function sublatticeMagnetization(spins::Array{Float64, 3}, n2Range::StepRange{Int64,Int64}, n3Range::StepRange{Int64,Int64})
  sublatticeM = [0., 0., 0.]
  for n2 in n2Range
    for n3 in n3Range
      sublatticeM += spins[n2, n3, :]
    end
  end
  sublatticeM
end

function MCRun(params::SystemParameters, bondMultipliers::Array{Float64, 3})
  system = Yb(params, bondMultipliers)  # initialize spins
  for sweep = 1:params.thermalizationSweeps  # thermalization sweeps
    MCSweep!(system, false)
  end

  # calculate the magnetization of each sublattice:
  evenRange = 2:2:params.L
  oddRange = 1:2:(params.L-1)
  system.MA = sublatticeMagnetization(system.spins, evenRange, evenRange)
  system.MB = sublatticeMagnetization(system.spins, oddRange, oddRange)
  system.MC = sublatticeMagnetization(system.spins, oddRange, evenRange)
  system.MD = sublatticeMagnetization(system.spins, evenRange, oddRange)

  # measurement sweeps:
  energyList = Float64[]
  psiList = Complex{Float64}[]  # measurements after each sweep
  for sweep = 1:params.equilibriumSweeps  # take equilibrium measurements
    MCSweep!(system, true)
    push!(energyList, system.energy)
    push!(psiList, (dot(system.MA, system.CA) + dot(system.MB, system.CB) + dot(system.MC, system.CC) + dot(system.MD, system.CD)) / system.N)
  end
  RawRunData(params.T, energyList / system.N, psiList)
end
