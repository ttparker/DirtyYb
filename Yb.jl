type Yb
  # system parameters:
  L::Int64
  T::Float64
  costheta::Float64  # these two only used in out-of-plane phase
  sintheta::Float64
  bondMultipliers::Array{Float64, 3} # (n_2, n_3, delta), where delta = {-a_1, a_2, a_3}
  # n.n. coupling matrices:
  Ja1::Array{Float64,2}
  Ja2::Array{Float64,2}
  Ja3::Array{Float64,2}

  spins::Array{Float64, 3}  # (n_2, n_3, s)
  energy::Float64

  # system-size-dependent quantities for making measurements:
  evenRange::StepRange{Int64,Int64}  # even integers from 2 to L
  oddRange::StepRange{Int64,Int64}  # odd integers from 1 to L-1
  N::Int64
end

function Yb(Jzz::Float64, Jpm::Float64, Jpmpm::Float64, Jzpm::Float64, bondMultipliers::Array{Float64, 3}, T::Float64, L::Int64)
  Ja1 = [2*(Jpm+Jpmpm) 0 0; 0 2*(Jpm-Jpmpm) Jzpm; 0 Jzpm Jzz]
  Ja2 = [2*Jpm-Jpmpm -sqrt(3)*Jpmpm -sqrt(3)/2*Jzpm; -sqrt(3)*Jpmpm 2*Jpm+Jpmpm -1/2*Jzpm; -sqrt(3)/2*Jzpm -1/2*Jzpm Jzz]
  Ja3 = [2*Jpm-Jpmpm sqrt(3)*Jpmpm sqrt(3)/2*Jzpm; sqrt(3)*Jpmpm 2*Jpm+Jpmpm -1/2*Jzpm; sqrt(3)/2*Jzpm -1/2*Jzpm Jzz]

  # initialize spins randomly:
  spins = Array{Float64}(L, L, 3)  # first index is n_2, second is n_3
  initialOrientations = rand(L, L, 2)
  for i in 1:L
    for j in 1:L
      spins[i, j, :] = randomSpin(initialOrientations[i, j, 1], initialOrientations[i, j, 2])
    end
  end
  # measure the energy of the initial spin configuration:
  energy = 0.
  for n2 in 1:L
    for n3 in 1:L
      nw = mod1(n2 + 1, L)
      sw = mod1(n3 + 1, L)
      energy += dot(spins[n2, n3, :], bondMultipliers[n2, n3, 1] * Ja1 * spins[nw, sw, :] + bondMultipliers[n2, n3, 2] * Ja2 * spins[nw, n3, :] + bondMultipliers[n2, n3, 3] * Ja3 * spins[n2, sw, :])  # only go "west" to avoid double-counting
    end
  end

  x = (Jpmpm + sqrt(Jzpm^2 + Jpmpm^2)) / Jzpm
  Yb(L, T, -1 / sqrt(1 + x^2), x / sqrt(1 + x^2), bondMultipliers, Ja1, Ja2, Ja3, spins, energy, 2:2:L, 1:2:(L-1), L^2)
end
