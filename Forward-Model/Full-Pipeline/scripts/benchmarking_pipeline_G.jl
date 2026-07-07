include("../../Neutrino-Flux/scripts/generate_neutrino_flux.jl")
using .Neutrino_Flux
include("../../Neutrino-Interactions/scripts/generate_neutrino_cross_sections.jl")
using .Neutrino_Cross_Sections
using Plots
using Interpolations
using LaTeXStrings

energies = logrange(1, 100, 100)
# read the neutrino-flux table nuflux.csv given in .../Neutrino-Flux/data
bin_centers, flux_νe_interp, flux_νμ_interp, flux_antiνe_interp, flux_antiνμ_interp, energies, flux_νe, flux_νμ, flux_antiνe, flux_antiνμ = Neutrino_Flux.read_neutrino_flux_table("nuflux", 100, 100, false)

# neutrino fluxes look good
# (comparing with Figure 3 in arXiv:1502.03916;
# there are slight differences, especially at low energies, because
# daemonflux is a bit different from honda)
p1 = plot(energies, (energies.^3) .* sum(flux_νe,     dims=2) ./ 99, label=L"$\nu_e$", 
          xaxis=:log, yaxis=:log, legend=:bottom,
          xlabel="E/GeV", ylabel=L"$\phi \times E_\nu^3$ /(m$^{-2}$s$^{-1}$sr$^{-1}$GeV$^{2}$")
plot!(p1, energies, (energies.^3) .* sum(flux_νμ,     dims=2) ./ 99, label=L"$\nu_\mu$")
plot!(p1, energies, (energies.^3) .* sum(flux_antiνe, dims=2) ./ 99, label=L"$\bar{\nu}_e$")
plot!(p1, energies, (energies.^3) .* sum(flux_antiνμ, dims=2) ./ 99, label=L"$\bar{\nu}_\mu$")

# plot the differences between the original and the interpolated fluxes
# to make sure that the interpolation is working correctly
# -> for the original flux we take the mean values
flux_νe_original     = 0.5 .* (sum(flux_νe,     dims=2)[1:end-1] .+ sum(flux_νe, dims=2)[2:end])
flux_νμ_original     = 0.5 .* (sum(flux_νμ,     dims=2)[1:end-1] .+ sum(flux_νμ, dims=2)[2:end])
flux_antiνe_original = 0.5 .* (sum(flux_antiνe, dims=2)[1:end-1] .+ sum(flux_antiνe, dims=2)[2:end])
flux_antiνμ_original = 0.5 .* (sum(flux_antiνμ, dims=2)[1:end-1] .+ sum(flux_antiνμ, dims=2)[2:end])
p2 = plot(bin_centers, 100 .* (sum(flux_νe_interp, dims=2) .- flux_νe_original) ./ flux_νe_original, 
          label=L"$\nu_e$", xlabel="E/GeV", ylabel="100 * (original - interpolation) / interpolation")
plot!(p2, bin_centers, 100 .* (sum(flux_νμ_interp,     dims=2) .- flux_νμ_original)     ./ flux_νμ_original,     label=L"$\nu_\mu$") 
plot!(p2, bin_centers, 100 .* (sum(flux_antiνe_interp, dims=2) .- flux_antiνe_original) ./ flux_antiνe_original, label=L"$\bar{\nu}_e$") 
plot!(p2, bin_centers, 100 .* (sum(flux_antiνμ_interp, dims=2) .- flux_antiνμ_original) ./ flux_antiνμ_original, label=L"$\bar{\nu}_\mu$") 





energies = logrange(0.1, 100, 100)
# read the neutrino osicllation probabilities calculated using Neutrhino





energies = logrange(0.1, 100, 100)
# read the neutrino cross sections from cross_sections.json in .../Neutrino-Interactions/data
cs_νe_CC, cs_νμ_CC, cs_ντ_CC, cs_antiνe_CC, cs_antiνμ_CC, cs_antiντ_CC = Neutrino_Cross_Sections.read_neutrino_cross_sections_info("cross_sections.json")
cs_νe_CC     = Neutrino_Cross_Sections.evaluate_cubicspline.(Ref(cs_νe_CC), log10.(energies))
cs_νμ_CC     = Neutrino_Cross_Sections.evaluate_cubicspline.(Ref(cs_νμ_CC), log10.(energies))
cs_ντ_CC     = Neutrino_Cross_Sections.evaluate_cubicspline.(Ref(cs_ντ_CC), log10.(energies))
cs_antiνe_CC = Neutrino_Cross_Sections.evaluate_cubicspline.(Ref(cs_antiνe_CC), log10.(energies))
cs_antiνμ_CC = Neutrino_Cross_Sections.evaluate_cubicspline.(Ref(cs_antiνμ_CC), log10.(energies))
cs_antiντ_CC = Neutrino_Cross_Sections.evaluate_cubicspline.(Ref(cs_antiντ_CC), log10.(energies))
# cross sections per nucleon look good 
# (comparing with EarthProbe/ExtModels/xsection/crossSection.root and
# with Figure 9 (top and bottom) in arXiv:1305.7513)
p3 = plot(energies, cs_νe_CC,     label=L"$\nu_e$ CC", xaxis=:log,
          xlabel="E/GeV", ylabel="(σ/E) 10^{-42} m^2/Gev")
plot!(p3, energies, cs_νμ_CC,     label=L"$\nu_\mu$ CC")
plot!(p3, energies, cs_ντ_CC,     label=L"$\nu_\tau$ CC")
plot!(p3, energies, cs_antiνe_CC, label=L"$\bar{\nu}_e$ CC")
plot!(p3, energies, cs_antiνμ_CC, label=L"$\bar{\nu}_\mu$ CC")
plot!(p3, energies, cs_antiντ_CC, label=L"$\bar{\nu}_\tau$ CC")


