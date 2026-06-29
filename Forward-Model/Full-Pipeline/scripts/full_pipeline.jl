include("../../Neutrino-Flux/scripts/create_neutrino_flux.jl")
using .Neutrino_Flux
include("../../Neutrino-Interactions/scripts/create_neutrino_cross_sections.jl")
using .Neutrino_Cross_Sections

energies = logrange(1, 100, 100)

# read the neutrino-flux table nuflux.csv given in .../Neutrino-Flux/data
flux_νe, flux_νμ, flux_antiνe, flux_antiνμ = read_neutrino_flux_table("nuflux", 100, 100, false)

# read the neutrino cross sections from cross_sections.json in .../Neutrino-Interactions/data
cs_νe_CC, cs_νμ_CC, cs_ντ_CC, cs_antiνe_CC, cs_antiνμ_CC, cs_antiντ_CC = read_neutrino_cross_sections_info("cross_sections.json")
cs_νe_CC = evaluate.(Ref(cs_νe_CC),  log10.(energies)) ./ energies)



