const dir = "JNeg/Delta0p4/L46/"
const filename = "Averaged.jld"

import JLD
using Plots
include("CustomTypes.jl")

function makePlot(plotyLabel::AbstractString, plotTitle::String, nRealizations::Int64, f::Function, system::SystemSummary, plotFilename::String)
  data = f(system)
  plot(xlabel = "T", ylabel = plotyLabel, title = plotTitle * " (" * string(nRealizations) * " realizations)\n" * (system.params.Jpmpm < 0 ? "In" : "Out-of") * "-plane phase, \$\\Delta = " * string(system.params.Delta) * "\$, \$N = " * string(system.params.L^2) * "\$\n")
  plot!(system.Ts, data, seriestype = :line, color = 1, label = "")
  plot!(system.Ts, data, seriestype = :scatter, color = 1, label = "")
  savefig("Analyzed/" * dir * "Plots/" * plotFilename)
end

function combinePlots(plotyLabel::AbstractString, plotTitle::String, nRealizations::Int64, fs::Vector{Function}, system::SystemSummary, labelList::Vector{String}, plotFilename::String)
  plot(xlabel = "T", ylabel = plotyLabel, title = plotTitle * " (" * string(nRealizations) * " realizations)\n" * (system.params.Jpmpm < 0 ? "In" : "Out-of") * "-plane phase, \$\\Delta = " * string(system.params.Delta) * "\$, \$N = " * string(system.params.L^2) * "\$\n")
  for i in 1:length(fs)
    data = fs[i](system)
    plot!(system.Ts, data, seriestype = :line, color = i, label = "")
    plot!(system.Ts, data, seriestype = :scatter, color = i, label = labelList[i])
  end
  savefig("Analyzed/" * dir * "Plots/" * plotFilename)
end

system = JLD.load("Analyzed/" * dir * filename, "system")::SystemSummary
nRealizations = JLD.load("Analyzed/" * dir * filename, "nRealizations")::Int64
mkpath("Analyzed/" * dir * "Plots/")

pyplot(markersize = 6)
Plots.scalefontsizes(1.5)
makePlot("E / N", "Energy density", nRealizations, x -> x.avgEnergies, system, "Energy.png")
makePlot("C", "Heat capacity", nRealizations, x -> x.energyVariances .* x.params.L.^2 ./ x.Ts.^2, system, "HeatCapacity.png")
makePlot("C", "Heat capacity", nRealizations, x -> x.energyVariances .* x.params.L.^2 ./ x.Ts.^2, system, "HeatCapacity.png")
makePlot("\$|\\langle \\psi(x) \\rangle|\$", "Order parameter", nRealizations, x -> x.absAvgPsis, system, "OrderParameter1.png")
makePlot("\$\\langle |\\psi(x)| \\rangle\$", "Order parameter", nRealizations, x -> x.avgAbsPsis, system, "OrderParameter2.png")
makePlot("\$|\\langle \\psi(x)^2 \\rangle|\$", "Order parameter", nRealizations, x -> x.absAvgPsi2s, system, "OrderParameter3.png")
makePlot("\$\\langle |\\psi(x)^2| \\rangle\$", "Order parameter", nRealizations, x -> x.avgAbsPsi2s, system, "OrderParameter4.png")
makePlot("\$|\\langle \\psi(x)^3 \\rangle|\$", "Order parameter", nRealizations, x -> x.absAvgPsi3s, system, "OrderParameter5.png")
makePlot("\$\\langle |\\psi(x)^3| \\rangle\$", "Order parameter", nRealizations, x -> x.avgAbsPsi3s, system, "OrderParameter6.png")

# Combine powers of order parameter:
combinePlots("\$\\left| \\left \\langle \\psi(x)^n \\right \\rangle \\right|\$", "Order parameters", nRealizations, [x -> x.absAvgPsis, x -> x.absAvgPsi2s, x -> x.absAvgPsi3s], system, ["\$\\left| \\left \\langle \\psi(x) \\right \\rangle \\right|\$", "\$\\left| \\left \\langle \\psi(x)^2 \\right \\rangle \\right|\$", "\$\\left| \\left \\langle \\psi(x)^3 \\right \\rangle \\right|\$"], "AbsAvg.png")

combinePlots("\$\\left \\langle \\left| \\psi(x)^n \\right| \\right \\rangle\$", "Order parameters", nRealizations, [x -> x.avgAbsPsis, x -> x.avgAbsPsi2s, x -> x.avgAbsPsi3s], system, ["\$\\left \\langle \\left| \\psi(x) \\right| \\right \\rangle\$", "\$\\left \\langle \\left| \\psi(x)^2 \\right| \\right \\rangle\$", "\$\\left \\langle \\left| \\psi(x)^3 \\right| \\right \\rangle\$"], "AvgAbs.png")

# Raised to powers to make comparable:
combinePlots("\$\\left| \\left \\langle \\psi(x)^n \\right \\rangle \\right|^{1/n}\$", "Order parameters", nRealizations, [x -> x.absAvgPsis, x -> sqrt.(x.absAvgPsi2s), x -> x.absAvgPsi3s.^(1/3)], system, ["\$\\left| \\left \\langle \\psi(x) \\right \\rangle \\right|\$", "\$\\left| \\left \\langle \\psi(x)^2 \\right \\rangle \\right|^{1/2}\$", "\$\\left| \\left \\langle \\psi(x)^3 \\right \\rangle \\right|^{1/3}\$"], "AbsAvgRaised.png")

combinePlots("\$\\left \\langle \\left| \\psi(x)^n \\right| \\right \\rangle^{1/n}\$", "Order parameters", nRealizations, [x -> x.avgAbsPsis, x -> sqrt.(x.avgAbsPsi2s), x -> x.avgAbsPsi3s.^(1/3)], system, ["\$\\left \\langle \\left| \\psi(x) \\right| \\right \\rangle\$", "\$\\left \\langle \\left| \\psi(x)^2 \\right| \\right \\rangle^{1/2}\$", "\$\\left \\langle \\left| \\psi(x)^3 \\right| \\right \\rangle^{1/3}\$"], "AvgAbsRaised.png")
