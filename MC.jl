function randomSpin(r1::Float64, r2::Float64)
  z = 2 * r1 - 1
  s = sqrt(1 - z^2)
  phi = 2 * pi * r2
  [s * cos(phi), s * sin(phi), z]
end

function calcPsi(system::Yb, spin::Vector{Float64}, n2::Int64, n3::Int64)
  n2Even = iseven(n2)
  n3Even = iseven(n3)
  if n2Even && n3Even
    C = system.CA
  elseif !n2Even && !n3Even
    C = system.CB
  elseif !n2Even && n3Even
    C = system.CC
  else
    C = system.CD
  end
  dot(spin, C)
end

function MCSweep!(system::Yb, measure::Bool)
  ns = rand(1:system.L, system.N, 2) # generating all rands at once saves RNG call overhead
  rands = rand(system.N, 3)  # (for z, for phi, for flipping prob)
  for step = 1:system.N
    n2, n3 = ns[step, :]  # pick a spin to flip
    nw = mod1(n2 + 1, system.L)
    se = mod1(n2 - 1, system.L)
    sw = mod1(n3 + 1, system.L)
    ne = mod1(n3 - 1, system.L)
    candidate = randomSpin(rands[step, 1], rands[step, 2])  # candidate new orientation
    DeltaE = dot(candidate - system.spins[n2, n3, :], system.Ja1 * (system.bondMultipliers[se, ne, 1] * system.spins[se, ne, :] + system.bondMultipliers[n2, n3, 1] * system.spins[nw, sw, :]) + system.Ja2 * (system.bondMultipliers[n2, n3, 2] * system.spins[nw, n3, :] + system.bondMultipliers[se, n3, 2] * system.spins[se, n3, :]) + system.Ja3 * (system.bondMultipliers[n2, n3, 3] * system.spins[n2, sw, :] + system.bondMultipliers[n2, ne, 3] * system.spins[n2, ne, :]))
    if DeltaE <= 0 || rands[step, 3] < exp(-DeltaE / system.T)
      if measure
        oldPsi = calcPsi(system, system.spins[n2, n3, :], n2, n3)
        newPsi = calcPsi(system, candidate, n2, n3)
        system.Psi += newPsi - oldPsi
        system.Psi2 += newPsi^2 - oldPsi^2
        system.Psi3 += newPsi^3 - oldPsi^3
      end
      system.spins[n2, n3, :] = candidate  # replace spin by candidate
      system.Energy += DeltaE
    end
  end
end

function MCRun(params::SystemParameters, bondMultipliers::Array{Float64, 3})
  system = Yb(params, bondMultipliers)  # initialize spins
  for sweep = 1:params.thermalizationSweeps  # thermalization sweeps
    MCSweep!(system, false)
  end

  # calculate the moments of psi:
  for n2 in 1:params.L
    for n3 in 1:params.L
      psi = calcPsi(system, system.spins[n2, n3, :], n2, n3)
      system.Psi += psi
      system.Psi2 += psi^2
      system.Psi3 += psi^3
    end
  end

  # measurement sweeps:
  EnergyList = Float64[]
  PsiList = Complex{Float64}[]  # measurements after each sweep
  Psi2List = Complex{Float64}[]
  Psi3List = Complex{Float64}[]

  for sweep = 1:params.equilibriumSweeps  # take equilibrium measurements
    MCSweep!(system, true)
    push!(EnergyList, system.Energy)
    push!(PsiList, system.Psi)
    push!(Psi2List, system.Psi2)
    push!(Psi3List, system.Psi3)
  end
  RawRunData(params.T, EnergyList / system.N, PsiList / system.N, Psi2List / system.N, Psi3List / system.N)
end
