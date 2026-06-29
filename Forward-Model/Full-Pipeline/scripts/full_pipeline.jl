include("../../Neutrino-Flux/scripts/create_neutrino_flux.jl")
using .Neutrino_Flux

# read the neutrino-flux table nuflux.csv given in .../Neutrino-Flux/data
nu_νe, nu_νμ, nu_antiνe, nu_antiνμ = read_neutrino_flux_table("nuflux", 100, 100, false)

