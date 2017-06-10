function MCRun(params::SystemParameters, bondMultipliers::Array{Float64, 3}, f::IOStream)
  system = Yb(params, bondMultipliers)  # initialize spins
  for sweep = 1:params.thermalizationSweeps  # wait for system to thermalize
    MCSweep!(system)
  end
  energyList = Float64[]
  psiList = Complex{Float64}[]  # measurements after each sweep
  for sweep = 1:params.equilibriumSweeps  # take equilibrium measurements
    MCSweep!(system)
    push!(energyList, system.energy)
    push!(psiList, measurePsi(system))
  end

  # write results to file:
  println(f, "T: ", params.T)
  println(f, "Energies:")
  writedlm(f, energyList' / system.N, ", ")
  println(f, "Psis:")
  writedlm(f, transpose(psiList), ", ")
  println(f)
end

function randomSpin(r1::Float64, r2::Float64)
  z = 2 * r1 - 1
  s = sqrt(1 - z^2)
  theta = 2 * pi * r2
  [s * cos(theta), s * sin(theta), z]
end

function MCSweep!(sys::Yb)
  ns = rand(1:sys.L, sys.N, 2) # generating all rands at once saves RNG call overhead
  rands = rand(sys.N, 3)
  for step = 1:sys.N
    n2, n3 = ns[step, :]  # pick a spin to flip
    nw = mod1(n2 + 1, sys.L)
    se = mod1(n2 - 1, sys.L)
    sw = mod1(n3 + 1, sys.L)
    ne = mod1(n3 - 1, sys.L)
    candidate = randomSpin(rands[step, 1], rands[step, 2])  # candidate new orientation
    DeltaE = dot(candidate - sys.spins[n2, n3, :], sys.Ja1 * (sys.bondMultipliers[se, ne, 1] * sys.spins[se, ne, :] + sys.bondMultipliers[n2, n3, 1] * sys.spins[nw, sw, :]) + sys.Ja2 * (sys.bondMultipliers[n2, n3, 2] * sys.spins[nw, n3, :] + sys.bondMultipliers[se, n3, 2] * sys.spins[se, n3, :]) + sys.Ja3 * (sys.bondMultipliers[n2, n3, 3] * sys.spins[n2, sw, :] + sys.bondMultipliers[n2, ne, 3] * sys.spins[n2, ne, :]))
    if DeltaE <= 0 || rands[step, 3] < exp(-DeltaE / sys.T)
      sys.spins[n2, n3, :] = candidate  # replace spin by candidate
      sys.energy += DeltaE
    end
  end
end

function sublatticeMagnetization(n2Range::StepRange{Int64,Int64}, n3Range::StepRange{Int64,Int64}, spins::Array{Float64, 3})
  sublatticeM = [0., 0., 0.]
  for n2 in n2Range
    for n3 in n3Range
      sublatticeM += spins[n2, n3, :]
    end
  end
  sublatticeM
end

function measurePsi(system::Yb)
  M1 = sublatticeMagnetization(system.evenRange, system.evenRange, system.spins)
  M2 = sublatticeMagnetization(system.oddRange, system.oddRange, system.spins)
  M3 = sublatticeMagnetization(system.oddRange, system.evenRange, system.spins)
  M4 = sublatticeMagnetization(system.evenRange, system.oddRange, system.spins)

  # in-plane phase (Jpmpm << 0):
  C1 = [3/2, 3/2*im, 0]
  C2 = [1/2, -3/2*im, 0]
  C3 = [-1-sqrt(3)/2*im, -sqrt(3)/2, 0]
  C4 = [-1+sqrt(3)/2*im, sqrt(3)/2, 0]

#=  # out-of-plane phase (Jpmpm >> 0) order parameter:
  C1 = [-3/2*im*system.sintheta, 3/2*system.sintheta, 0]
  C2 = [3/2*im*system.sintheta, 1/2*system.sintheta, 2*system.costheta]
  C3 = [sqrt(3)/2*system.sintheta, (-1-sqrt(3)/2*im)*system.sintheta, (-1+sqrt(3)*im)*system.costheta]
  C4 = [-sqrt(3)/2*system.sintheta, (-1+sqrt(3)/2*im)*system.sintheta, (-1-sqrt(3)*im)*system.costheta] =#

  (dot(M1, C1) + dot(M2, C2) + dot(M3, C3) + dot(M4, C4)) / system.N
end
