if isempty(ARGS)
  const L = 10
  const Delta = 0.2
else # args from command line override those values
  const L = parse(Int64, ARGS[1])
  const Delta = parse(Float64, ARGS[2])
end

import JLD

JLD.save("Realization.jld", "L", L, "Delta", Delta, "multipliers", fill(1 - Delta, L, L, 3) + 2 * Delta * rand(L, L, 3))
