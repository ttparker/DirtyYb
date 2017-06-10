const L = 10
const Delta = 0.2

import JLD

JLD.save("Realization.jld", "L", L, "Delta", Delta, "multipliers", fill(1 - Delta, L, L, 3) + 2 * Delta * rand(L, L, 3))
