include("../../Neutrino-Flux/scripts/create_neutrino_flux.jl")
using .Neutrino_Flux
include("../../Neutrino-Interactions/scripts/create_neutrino_cross_sections.jl")
using .Neutrino_Cross_Sections
using Plots
using Interpolations

energies = logrange(1, 100, 100)
# read the neutrino-flux table nuflux.csv given in .../Neutrino-Flux/data
flux_νe, flux_νμ, flux_antiνe, flux_antiνμ = read_neutrino_flux_table("nuflux", 100, 100, false)
# neutrino fluxes look good
# (comparing with Figure 3 in 2arXiv:1502.03916,
# there are slight differences, especially at low energies, because
# daemonflux is a bit different from honda)
p1 = plot(energies, (energies.^3) .* sum(flux_νe,     dims=2) ./ 99, label="νe", 
          xaxis=:log, yaxis=:log, legend=:bottom,
          xlabel="E/GeV", ylabel="(σ/E) 10^{-42} m^2/Gev")
plot!(p1, energies, (energies.^3) .* sum(flux_νμ,     dims=2) ./ 99, label="νμ")
plot!(p1, energies, (energies.^3) .* sum(flux_antiνe, dims=2) ./ 99, label="antiνe")
plot!(p1, energies, (energies.^3) .* sum(flux_antiνμ, dims=2) ./ 99, label="antiνμ")



energies = logrange(0.1, 100, 100)
# read the neutrino cross sections from cross_sections.json in .../Neutrino-Interactions/data
cs_νe_CC, cs_νμ_CC, cs_ντ_CC, cs_antiνe_CC, cs_antiνμ_CC, cs_antiντ_CC = read_neutrino_cross_sections_info("cross_sections.json")
cs_νe_CC     = evaluate_cubicspline.(Ref(cs_νe_CC), log10.(energies))
cs_νμ_CC     = evaluate_cubicspline.(Ref(cs_νμ_CC), log10.(energies))
cs_ντ_CC     = evaluate_cubicspline.(Ref(cs_ντ_CC), log10.(energies))
cs_antiνe_CC = evaluate_cubicspline.(Ref(cs_antiνe_CC), log10.(energies))
cs_antiνμ_CC = evaluate_cubicspline.(Ref(cs_antiνμ_CC), log10.(energies))
cs_antiντ_CC = evaluate_cubicspline.(Ref(cs_antiντ_CC), log10.(energies))
# cross sections per nucleon look good 
# (comparing with EarthProbe/ExtModels/xsection/crossSection.root and
# with Figure 9 (top and bottom) in arXiv:1305.7513)
p2 = plot(energies, cs_νe_CC,     label="νe CC", xaxis=:log,
          xlabel="E/GeV", ylabel="(σ/E) 10^{-42} m^2/Gev")
plot!(p2, energies, cs_νμ_CC,     label="νμ CC")
plot!(p2, energies, cs_ντ_CC,     label="ντ CC")
plot!(p2, energies, cs_antiνe_CC, label="anti-νe CC")
plot!(p2, energies, cs_antiνμ_CC, label="anti-νμ CC")
plot!(p2, energies, cs_antiντ_CC, label="anti-ντ CC")