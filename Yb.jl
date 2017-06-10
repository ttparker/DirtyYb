type Yb
  # system parameters:
  L::Int64
  T::Float64
  costheta::Float64  # these two only used in out-of-plane phase
  sintheta::Float64
  bondMultipliers::Array{Float64, 3} # (n_2, n_3, delta), where delta is in {-a_1, a_2, a_3}
  # n.n. coupling matrices:
  Ja1::Matrix{Float64}
  Ja2::Matrix{Float64}
  Ja3::Matrix{Float64}

  spins::Array{Float64, 3}  # (n_2, n_3, s)
  energy::Float64

  # system-size-dependent quantities for making measurements:
  evenRange::StepRange{Int64,Int64}  # even integers from 2 to L
  oddRange::StepRange{Int64,Int64}  # odd integers from 1 to L-1
  N::Int64
end

function Yb(params::SystemParameters, bondMultipliers::Array{Float64, 3})
  Ja1 = [2*(params.Jpm+params.Jpmpm) 0 0; 0 2*(params.Jpm-params.Jpmpm) params.Jzpm; 0 params.Jzpm params.Jzz]
  Ja2 = [2*params.Jpm-params.Jpmpm -sqrt(3)*params.Jpmpm -sqrt(3)/2*params.Jzpm; -sqrt(3)*params.Jpmpm 2*params.Jpm+params.Jpmpm -1/2*params.Jzpm; -sqrt(3)/2*params.Jzpm -1/2*params.Jzpm params.Jzz]
  Ja3 = [2*params.Jpm-params.Jpmpm sqrt(3)*params.Jpmpm sqrt(3)/2*params.Jzpm; sqrt(3)*params.Jpmpm 2*params.Jpm+params.Jpmpm -1/2*params.Jzpm; sqrt(3)/2*params.Jzpm -1/2*params.Jzpm params.Jzz]

  # initialize spins randomly:
  spins = Array{Float64, 3}(params.L, params.L, 3)
  initialOrientations = rand(params.L, params.L, 2)
  for i in 1:params.L
    for j in 1:params.L
      spins[i, j, :] = randomSpin(initialOrientations[i, j, 1], initialOrientations[i, j, 2])
    end
  end
  # measure the energy of the initial spin configuration:
  energy = 0.
  for n2 in 1:params.L
    for n3 in 1:params.L
      nw = mod1(n2 + 1, params.L)
      sw = mod1(n3 + 1, params.L)
      energy += dot(spins[n2, n3, :], bondMultipliers[n2, n3, 1] * Ja1 * spins[nw, sw, :] + bondMultipliers[n2, n3, 2] * Ja2 * spins[nw, n3, :] + bondMultipliers[n2, n3, 3] * Ja3 * spins[n2, sw, :])  # only go "west" to avoid double-counting
    end
  end

  x = (params.Jpmpm + sqrt(params.Jzpm^2 + params.Jpmpm^2)) / params.Jzpm
  Yb(params.L, params.T, -1 / sqrt(1 + x^2), x / sqrt(1 + x^2), bondMultipliers, Ja1, Ja2, Ja3, spins, energy, 2:2:params.L, 1:2:(params.L-1), params.L^2)
end
