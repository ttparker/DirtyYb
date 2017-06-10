const filenames = ["L46", "L64"]  # Names of analyzed data files to combine into plots

import JLD
include("CustomTypes.jl")

const nSystems = length(filenames)
function combineSystems(plotyLabel::AbstractString, plotTitle::String, f::Function, systems::Vector{SystemSummary}, plotFilename::String)
  plot(xlabel = "T", ylabel = plotyLabel, title = plotTitle)
  for i in 1:nSystems
    data = f(systems[i])
    plot!(systems[i].Ts, data, seriestype = :line, color = i, label = "")
    plot!(systems[i].Ts, data, seriestype = :scatter, color = i, label = systems[i].legendLabel)
  end
  savefig("Analyzed/" * plotFilename * ".png")
end

# Read in data from Analyzed folder:
if(!isdir("Analyzed/"))
  println("No AnalyzedData directory found - run DataAnalysis.jl first.")
  quit()
end
systems = Vector{SystemSummary}()
for filename in filenames
  push!(systems, JLD.load("Analyzed/" * filename * ".jld", filename))
end

# Make plots:
using Plots

pyplot(markersize = 6)
Plots.scalefontsizes(1.75)
combineSystems("E/N", "Energy density", x -> x.avgEnergies, systems, "Energy")
combineSystems("C", "Heat capacity", x -> x.energyVariances .* x.params.L.^2 ./ x.Ts.^2, systems, "HeatCapacity")
combineSystems("\$|\\langle \\psi \\rangle|\$", "Order parameter", x -> x.absAvgPsis, systems, "OrderParameter")
combineSystems("\$Q_2\$", "Binder ratio \$Q_2\$", x -> x.avgAbsPsi2s ./ x.avgAbsPsis.^2, systems, "BinderRatio")

# Break down various order parameters for each system:
for system in systems
  plot(system.Ts, [system.absAvgPsis, system.avgAbsPsis, sqrt(system.absAvgPsi2s), sqrt(system.avgAbsPsi2s)], seriestype = :line, label = "", xlabel = "T", title = "Order parameter for " * system.legendLabel)
  plot!(system.Ts, system.absAvgPsis, seriestype = :scatter, color = 1, markershape = :circle, label = "\$|\\langle \\psi \\rangle|\$")
  plot!(system.Ts, system.avgAbsPsis, seriestype = :scatter, color = 2, markershape = :rect, label = "\$\\langle |\\psi| \\rangle\$")
  plot!(system.Ts, sqrt(system.absAvgPsi2s), seriestype = :scatter, color = 3, markershape = :diamond, label = "\$\\sqrt{|\\langle \\psi^2 \\rangle|}\$")
  plot!(system.Ts, sqrt(system.avgAbsPsi2s), seriestype = :scatter, color = 4, markershape = :star5, label = "\$\\sqrt{\\langle |\\psi^2| \\rangle}\$")
  savefig("Analyzed/L" * string(system.params.L) * "OrderParameter.png")
end
