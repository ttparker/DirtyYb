const dir = "Delta0p4/L46/"
const realizations = ["Real1", "Real2", "Real3", "Real4", "Real5", "Real6"]

import JLD
using Plots
include("CustomTypes.jl")
pyplot(markersize = 6)
Plots.scalefontsizes(1.75)

function drawPlot(system::SystemSummary, f::Function, plotyLabel::AbstractString, plotTitle::String, plotFilename::String)
  plot(xlabel = "T", ylabel = plotyLabel, title = plotTitle)
  data = f(system)
  plot!(system.Ts, data, seriestype = :line, color = 1, label = "")
  plot!(system.Ts, data, seriestype = :scatter, color = 1, label = system.legendLabel)
  savefig("Analyzed/" * dir * plotFilename * ".png")
end

for realization in realizations
  # Read in data from Analyzed folder:
  if(!isdir("Analyzed/" * dir))
    println("No AnalyzedData directory found - run DataAnalysis.jl first.")
    quit()
  end
  system = JLD.load("Analyzed/" * dir * realization * "/" * realization * ".jld", "system")

  # Make plots:
  drawPlot(system, x -> x.avgEnergies, "E/N", "Energy density", realization * "/Energy")
  drawPlot(system, x -> x.energyVariances .* x.params.L.^2 ./ x.Ts.^2, "C", "Heat capacity", realization * "/HeatCapacity")

  # Break down various order parameters:
  plot(system.Ts, [system.absAvgPsis, system.avgAbsPsis, sqrt(system.absAvgPsi2s), sqrt(system.avgAbsPsi2s)], seriestype = :line, label = "", xlabel = "T", title = "Order parameter for " * system.legendLabel)
  plot!(system.Ts, system.absAvgPsis, seriestype = :scatter, color = 1, markershape = :circle, label = "\$|\\langle \\psi \\rangle|\$")
  plot!(system.Ts, system.avgAbsPsis, seriestype = :scatter, color = 2, markershape = :rect, label = "\$\\langle |\\psi| \\rangle\$")
  plot!(system.Ts, sqrt(system.absAvgPsi2s), seriestype = :scatter, color = 3, markershape = :diamond, label = "\$\\sqrt{|\\langle \\psi^2 \\rangle|}\$")
  plot!(system.Ts, sqrt(system.avgAbsPsi2s), seriestype = :scatter, color = 4, markershape = :star5, label = "\$\\sqrt{\\langle |\\psi^2| \\rangle}\$")
  savefig("Analyzed/" * dir * realization * "/OrderParameter.png")
end
