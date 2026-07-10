using  Pkg

#%% Relevant Directories
const FLEXOPT_DIR = "/Users/igoos/Desktop/projects/flexOPT" # TO-DO: adapt this path to FlexOPT

const FORWARDMODEL_DIR = "/Users/igoos/Desktop/projects/Joint-Inversion/Forward-Model" # TO-DO> adapt this path
#cd(FLEXOPT_DIR)

#%% 1. Activate the flexOPT environment
Pkg.activate(FLEXOPT_DIR)

#%% 2. Standard packages for plotting
using CairoMakie
using DIVAnd
using Interpolations
using Colors
using DelimitedFiles
using JSON
#using FilePaths
using FilePathsBase

#%% 3. Pipeline modules
include(joinpath(FORWARDMODEL_DIR, "Neutrino-Flux","scripts", "create_neutrino_flux.jl"))
using .Neutrino_Flux

#include("../../Neutrino-Interactions/scripts/create_neutrino_cross_sections.jl")
include(joinpath(FORWARDMODEL_DIR, "Neutrino-Interactions","scripts", "create_neutrino_cross_sections.jl"))
using .Neutrino_Cross_Sections

#include("../../Neutrino-Interactions/scripts/create_neutrino_cross_sections.jl")
ParamFile = "../config/testparam.csv" 
include(joinpath(FORWARDMODEL_DIR, "Neutrino-Probabilities","scripts", "generate_oscillation_probabilities.jl"))
using .Neutrino_Oscillation



#%%

# NOTES FROM YAEL: Hey Isa here some notes: 
#=
    1) I wanted to calculate the bin centers for bins evenly space in log10(E). 
    then i realized when you generate the flux you dont evaluate the flux at the bin centers. but
    rather you generate the flux tables using the points at which the flux is evaluated. 
    I will leave my first draft of the bin centers and widths below, but i will keep evaluated as you did in the original
    version of this file.

    2) using plot is extremely conclicting with the other packages, so i removed it an will use CairoMakie instead for all plotting.
    porpuses, including the plotting of the fluxes, cross sections and oscillation probabilities (HeatMaps).

 =#

# BIN CENTER FOR LOG10(E) EVENLY SPACED ENERGY BINS: ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Upper edges of the energy bins in GeV, In this case they are evenly spaced in log10(E) from 1 GeV to 100 GeV
E_uedges = logrange(1, 100, 100) # WARNING, FIRST ELEMENT IS THE LOWER EDGE OF THE FIRST BIN, NOT THE UPPER EDGE
energies = 0.5 .* (E_uedges[1:end-1] .+ E_uedges[2:end]) #Energy bin centers in GeV
ΔE = diff(energies)

fig = Figure(size = (800, 400))
ax  = Axis(fig[1, 1], 
           xlabel = "Energy (GeV)",
           title  = "Energy bin edges, centers and widths")

# lower and upper edges as vertical lines
vlines!(ax, E_uedges, color = :black, linestyle = :dash, label = "edges")

# centers as points
scatter!(ax, energies, fill(1.0, length(energies)),
         color = :red, markersize = 12, label = "centers")

# bin widths as horizontal spans
for i in 1:length(energies)
    vspan!(ax, E_uedges[i], E_uedges[i+1],
           color = (:blue, 0.1 + 0.1 * (i % 2)))  # alternating shading
end

axislegend(ax, position = :lt)
ylims!(ax, 0, 2)        # y axis is meaningless here, just for visual space
hideydecorations!(ax)   # hide y axis since it carries no information
display(fig)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#%% Pipeline variables ------------------------------------------------------------------------------------------------------------------------------

# From original benchmarking_pipeline_G.jl file
energies = logrange(1, 100, 100)
ΔE = diff(energies) # We need 100 widths, so for now i will include tthe last width twice, but we should fix this in the future
#energies = 10 .^ range(0, stop=2, length=100)
n_vectors  = 100       # WARNING HARDCODED IN Neutrino-Probabilities/scripts/generate_oscillation_probabilities.jl -> Must be the same as in generate_oscillation_probabilities.jl file
cos_θ    = range(-1, 0, length = n_vectors)
Δcosθ = diff(cos_θ)[1] 
#remember: fluxes are given in (m^2 sec sr GeV)^-1
time  = 365.25 * 24.0 * 60.0 * 60.0 #year in seconds -> I want 1Mton*year exposure to compare with Simon
p_m   = 1.67262 * 1000000.0 #+42-36 #proton mass in Mton
det_m = 1.0 #in Mton 


ξ_Det = time * det_m # Detector exposure in sec * Mton 


#%% Check: Plot Neutrino Flux------------------------------------------------------------------------------------------------------------------------------

energies = logrange(1, 100, 100)
ΔE = diff(energies)
push!(ΔE, energies[end] * (energies[end] / energies[end-1] - 1))
# read the neutrino-flux table nuflux.csv given in .../Neutrino-Flux/data
flux_νe, flux_νμ, flux_antiνe, flux_antiνμ = Neutrino_Flux.read_neutrino_flux_table("nuflux", 100, 100, false)
# neutrino fluxes look good
# (comparing with Figure 3 in 2arXiv:1502.03916;
# there are slight differences, especially at low energies, because
# daemonflux is a bit different from honda)
# remove: using Plots
# replace p1 = plot(...) with:
fig2, ax2, l2 = lines(log10.(energies), log10.((energies.^3) .* sum(flux_νe, dims=2)[:] ./ 99), 
                    label="νe", axis=(xlabel="E/GeV", ylabel="E³Φ"))
lines!(ax2, log10.(energies), log10.((energies.^3) .* sum(flux_νμ, dims=2)[:] ./ 99), label="νμ")

lines!(ax2, log10.(energies), log10.((energies.^3) .* sum(flux_antiνe, dims=2)[:] ./ 99), label="antiνe")
lines!(ax2, log10.(energies), log10.((energies.^3) .* sum(flux_antiνμ, dims=2)[:] ./ 99), label="antiνμ")

axislegend(ax2, position=:lb)
display(fig2)

#%% Neutrino Cross Sections check------------------------------------------------------------------------------------------------------------------------------

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
fig3, ax3, _ = lines(energies, cs_νe_CC ./ energies,
                     label = "νe CC",
                     axis  = (xscale = log10,
                              xlabel = "E/GeV",
                              ylabel = "(σ/E) 10⁻⁴² m²/GeV"))
lines!(ax3, energies, cs_νμ_CC     ./ energies, label = "νμ CC")
lines!(ax3, energies, cs_ντ_CC     ./ energies, label = "ντ CC")
lines!(ax3, energies, cs_antiνe_CC ./ energies, label = "anti-νe CC")
lines!(ax3, energies, cs_antiνμ_CC ./ energies, label = "anti-νμ CC")
lines!(ax3, energies, cs_antiντ_CC ./ energies, label = "anti-ντ CC")
axislegend(ax3, position = :rb)
display(fig3)


#%%  CHECK Neutrino Oscillation Probabilities plot------------------------------------------------------------------------------------------------------------------------------
energies = collect(logrange(1, 100, 100)) # WARNING: Logrange need to be "collected" to be used in the linkWithNeurthinoPREM_YD function
probs2, probs2anti = Neutrino_Oscillation.Generate_Oscillation_Probabilities(energies)

# An example of how to plot the oscillation probability P(νe → νμ) 
Pemu =Matrix(probs2[:,:,1,2])

fig = Figure()
ax  = Axis(fig[1, 1], aspect = 1, xscale = log10,
           xlabel = "Energy (GeV)", ylabel = "cos(θ)")
hm  = CairoMakie.heatmap!(ax, energies, collect(cos_θ), Pemu,
                           colormap = :inferno)
xlims!(ax, extrema(energies))
ylims!(ax, -1, 0)
Colorbar(fig[:, 2], hm, label = "Probability")
display(fig)

#%% Full interactive pipeline check------------------------------------------------------------------------------------------------------------------------------

int_evts_νe     =  (ξ_Det/p_m) * 2.0 * π * Δcosθ .* (cs_νe_CC) .* (flux_νe     .* probs2[:, :, 1, 1]     .+ flux_νμ     .* probs2[:, :, 2, 1])     .* ΔE
int_evts_νμ     =  (ξ_Det/p_m) * 2.0 * π * Δcosθ .* (cs_νμ_CC) .* (flux_νe     .* probs2[:, :, 1, 2]     .+ flux_νμ     .* probs2[:, :, 2, 2])     .* ΔE
int_evts_ντ     =  (ξ_Det/p_m) * 2.0 * π * Δcosθ .* (cs_ντ_CC) .* (flux_νe     .* probs2[:, :, 1, 3]     .+ flux_νμ     .* probs2[:, :, 2, 3])     .* ΔE
int_evts_antiνe =  (ξ_Det/p_m) * 2.0 * π * Δcosθ .* (cs_antiνe_CC) .* (flux_antiνe .* probs2anti[:, :, 1, 1] .+ flux_antiνμ .* probs2anti[:, :, 2, 1]) .* ΔE
int_evts_antiνμ =  (ξ_Det/p_m) * 2.0 * π * Δcosθ .* (cs_antiνμ_CC) .* (flux_antiνe .* probs2anti[:, :, 1, 2] .+ flux_antiνμ .* probs2anti[:, :, 2, 2]) .* ΔE
int_evts_antiντ =  (ξ_Det/p_m) * 2.0 * π * Δcosθ .* (cs_antiντ_CC) .* (flux_antiνe .* probs2anti[:, :, 1, 3] .+ flux_antiνμ .* probs2anti[:, :, 2, 3]) .* ΔE


# data and labels
channels = [
    (int_evts_νe,     "νe CC"),
    (int_evts_νμ,     "νμ CC"),
    (int_evts_ντ,     "ντ CC"),
    (int_evts_antiνe, "anti-νe CC"),
    (int_evts_antiνμ, "anti-νμ CC"),
    (int_evts_antiντ, "anti-ντ CC"),
]

fig5 = Figure(size = (1200, 700))

for (idx, (data, title)) in enumerate(channels)
    row = (idx - 1) ÷ 3 + 1   # rows 1-2
    col = (idx - 1) % 3 + 1   # cols 1-3

    ax = Axis(fig5[row, col],
              title  = title,
              xscale = log10,
              xlabel = "Energy (GeV)",
              ylabel = "cos(θ)")

    hm = CairoMakie.heatmap!(ax, energies, collect(cos_θ), parent(data),
                              colormap = :inferno)

    xlims!(ax, extrema(energies))
    ylims!(ax, -1, 0)
    Colorbar(fig5[row, col + 3], hm)   # colorbar to the right of each column... 
end                                     # ...or use a shared one (see below)

Label(fig5[0, :], "Interacting Events per bin", fontsize = 18)
display(fig5)

#%%

function diff(probs, probs2)
    cos_θ = range(-1, 0, length = n_vectors)
    fig3 = Figure()

    probdiff = probs-probs2
    matprobdiff = parent(probdiff)
    ax3 = Axis(fig3[1,1], aspect = 1, xscale=log10, xlabel="Energy (GeV)", ylabel="cos(θ)")
    hm3=heatmap!(ax3, energies, cos_θ, matprobdiff, colormap=cgrad(:inferno))
    Colorbar(fig3[:,2], hm3, label="Probability")
    display(fig3)
end
