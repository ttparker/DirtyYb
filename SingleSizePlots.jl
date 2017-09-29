const filename = "JNeg/Delta0p4/L46/Averaged.jld"

import JLD
using Plots
include("CustomTypes.jl")

function makePlot(plotyLabel::AbstractString, plotTitle::String, nRealizations::Int64, f::Function, system::SystemSummary, plotFilename::String)
  data = f(system)
  plot(xlabel = "T", ylabel = plotyLabel, title = plotTitle * " (" * string(nRealizations) * " realizations)\n" * (system.params.Jpmpm < 0 ? "In" : "Out-of") * "-plane phase, \$\\Delta = " * string(system.params.Delta) * "\$, \$N = " * string(system.params.L^2) * "\$\n")
  plot!(system.Ts, data, seriestype = :line, color = 1, label = "")
  plot!(system.Ts, data, seriestype = :scatter, color = 1, label = "")
  savefig("SingleSizePlots/" * plotFilename * ".png")
end

mkpath("SingleSizePlots/")
system = JLD.load("Analyzed/" * filename, "system")::SystemSummary
nRealizations = JLD.load("Analyzed/" * filename, "nRealizations")::Int64

pyplot(markersize = 6)
Plots.scalefontsizes(1.5)
makePlot("E / N", "Energy density", nRealizations, x -> x.avgEnergies, system, "Energy")
makePlot("C", "Heat capacity", nRealizations, x -> x.energyVariances .* x.params.L.^2 ./ x.Ts.^2, system, "HeatCapacity")
makePlot("C", "Heat capacity", nRealizations, x -> x.energyVariances .* x.params.L.^2 ./ x.Ts.^2, system, "HeatCapacity")
makePlot("\$|\\langle \\psi(x) \\rangle|\$", "Order parameter", nRealizations, x -> x.absAvgPsis, system, "OrderParameter1")
makePlot("\$\\langle |\\psi(x)| \\rangle\$", "Order parameter", nRealizations, x -> x.avgAbsPsis, system, "OrderParameter2")
makePlot("\$|\\langle \\psi(x)^2 \\rangle|\$", "Order parameter", nRealizations, x -> x.absAvgPsi2s, system, "OrderParameter3")
makePlot("\$\\langle |\\psi(x)^2| \\rangle\$", "Order parameter", nRealizations, x -> x.avgAbsPsi2s, system, "OrderParameter4")
makePlot("\$|\\langle \\psi(x)^3 \\rangle|\$", "Order parameter", nRealizations, x -> x.absAvgPsi3s, system, "OrderParameter5")
makePlot("\$\\langle |\\psi(x)^3| \\rangle\$", "Order parameter", nRealizations, x -> x.avgAbsPsi3s, system, "OrderParameter6")

# Combine powers of order parameter:

plot(xlabel = "T", ylabel = "\$|\\langle \\psi(x)^n \\rangle |\$", title = "Order parameters (" * string(nRealizations) * " realizations)\n" * (system.params.Jpmpm < 0 ? "In" : "Out-of") * "-plane phase, \$\\Delta = " * string(system.params.Delta) * "\$, \$N = " * string(system.params.L^2) * "\$\n")
plot!(system.Ts, system.absAvgPsis, seriestype = :line, color = 1, label = "")
plot!(system.Ts, system.absAvgPsis, seriestype = :scatter, color = 1, label = "\$|\\langle \\psi(x) \\rangle|\$")
plot!(system.Ts, system.absAvgPsi2s, seriestype = :line, color = 2, label = "")
plot!(system.Ts, system.absAvgPsi2s, seriestype = :scatter, color = 2, label = "\$|\\langle \\psi(x)^2 \\rangle|\$")
plot!(system.Ts, system.absAvgPsi3s, seriestype = :line, color = 3, label = "")
plot!(system.Ts, system.absAvgPsi3s, seriestype = :scatter, color = 3, label = "\$|\\langle \\psi(x)^3 \\rangle|\$")
savefig("SingleSizePlots/AbsAvg.png")

plot(xlabel = "T", ylabel = "\$\\langle |\\psi(x)^n| \\rangle \$", title = "Order parameters (" * string(nRealizations) * " realizations)\n" * (system.params.Jpmpm < 0 ? "In" : "Out-of") * "-plane phase, \$\\Delta = " * string(system.params.Delta) * "\$, \$N = " * string(system.params.L^2) * "\$\n")
plot!(system.Ts, system.avgAbsPsis, seriestype = :line, color = 1, label = "")
plot!(system.Ts, system.avgAbsPsis, seriestype = :scatter, color = 1, label = "\$\\langle |\\psi(x)| \\rangle\$")
plot!(system.Ts, system.avgAbsPsi2s, seriestype = :line, color = 2, label = "")
plot!(system.Ts, system.avgAbsPsi2s, seriestype = :scatter, color = 2, label = "\$\\langle |\\psi(x)^2| \\rangle\$")
plot!(system.Ts, system.avgAbsPsi3s, seriestype = :line, color = 3, label = "")
plot!(system.Ts, system.avgAbsPsi3s, seriestype = :scatter, color = 3, label = "\$\\langle |\\psi(x)^3| \\rangle\$")
savefig("SingleSizePlots/AvgAbs.png")
