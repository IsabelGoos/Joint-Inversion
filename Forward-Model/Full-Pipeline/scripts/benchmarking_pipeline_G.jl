#using  Pkg
#Pkg.activate(FLEXOPT_DIR)
#using DotEnv
#DotEnv.load!()
#const FLEXOPT_DIR = ENV["FLEXOPT_DIR"]
#using Pkg
#Pkg.activate(FLEXOPT_DIR)
#const FLEXOPT_DIR = "/Users/igoos/Desktop/projects/flexOPT" # TO-DO: adapt this path to FlexOPT

include("../../Neutrino-Flux/scripts/generate_neutrino_flux.jl")
using .Neutrino_Flux
ParamFile = "../config/testparam.csv" # TO-DO: adapt this path, in case it changes 
include("../../Neutrino-Oscillations/scripts/generate_oscillation_probabilities.jl")
using .Neutrino_Oscillations
include("../../Neutrino-Interactions/scripts/generate_neutrino_cross_sections.jl")
using .Neutrino_Cross_Sections
using Plots
using Interpolations
using LaTeXStrings #TO-DO: Add to manifesto &project

energies = logrange(1, 100, 100)
energies = logrange(1, 40, 100)
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

# This is what we need to compute the interacting events!:
# TO-DO: put this computation in the function, no user should be doing this XD
flux_νe = reverse(flux_νe', dims=1)
flux_νμ = reverse(flux_νμ', dims=1)
flux_antiνe = reverse(flux_antiνe', dims=1)
flux_antiνμ = reverse(flux_antiνμ', dims=1)




egrid = logrange(1, 40, 1000)
cosθgrid = range(-1, 0, 100)
nuflux_params_dict = Dict(
    "model"                => "Daemonflux",
    "bin_centers_arrays"   => (egrid, cosθgrid),
    "flux_mode"            => "conventional",
    "return_uncertainties" => false
)
output_dir = Neutrino_Flux.produce_neutrino_flux(nuflux_params_dict)
NuMu_flux = output_dir["NuMu_flux"]
p = heatmap(egrid, cosθgrid, log10.(NuMu_flux), xaxis=:log)




using Downloads
using CodecZlib
io_buffer = IOBuffer()
url = "http://www-rccn.icrr.u-tokyo.ac.jp/mhonda/public/nflx2014/kam-1202-20-01-solmin.d.gz"
Downloads.download(url, io_buffer)
seekstart(io_buffer)
decompressed_stream = GzipDecompressorStream(io_buffer)
lines = readlines(decompressed_stream)
close(decompressed_stream)
numeric_lines = filter(line -> occursin(r"^\s*[0-9.-]", line), lines)
raw_data = map(numeric_lines) do line
    parse.(Float64, split(line))
end
matrix_2d = hcat(raw_data...)
n_E = 101
n_cosθ = 20
n_azimuth = 12 

flux_grid = reshape(matrix_2d, 5, n_E, n_cosθ)
energies   = flux_grid[1, :, 1]
numu_flux   = flux_grid[2, :, :]
anumu_flux  = flux_grid[3, :, :]
nue_flux    = flux_grid[4, :, :]
anue_flux   = flux_grid[5, :, :]
cosθgrid = range(-1, 1, 20)
p = heatmap(energies, cosθgrid, log10.(numu_flux'), xscale=:log10)
size(energies)
size(cosθgrid)
size(numu_flux)
for i = 1:length(energies)
    println(energies[i])  # Changed 1 to i
end

nuflux_params = Dict{String, Any}(
    "model"                  => :Honda,
    "location_hf"            => "kamioka",      
    "season_hf"              => "dec_feb",
    "angles_separation_hf"   => "averaged_ϕ",       
    "mountain_overburden_hf" => "without",
    "solar_activity_hf"      => "minimum"   
)
nufluxresult = produce_neutrino_flux(nuflux_params)

energies = nufluxresult["Energies"]
cosθgrid = nufluxresult["cosθs"]
numu_flux = nufluxresult["NuMu_flux"]
nue_flux = nufluxresult["Nue_flux"]
anumu_flux = nufluxresult["AntiNuMu_flux"]

p = heatmap(energies, cosθgrid, log10.(numu_flux), xscale=:log10)

my_flux = Neutrino_Flux.DEFAULT_NUFLUX[]
def_energies = my_flux["Energies"]
def_cosθgrid = my_flux["cosθs"]
def_numu_flux = my_flux["NuMu_flux"]
def_nue_flux = my_flux["Nue_flux"]
def_anumu_flux = my_flux["AntiNuMu_flux"]

p = heatmap(def_energies, def_cosθgrid, log10.(def_numu_flux), xscale=:log10)




# neutrino oscillation probabilities look good
# TO-DO: Paνe2aντ deserves some attention → do a difference with OscProb!
minX, maxX, nX = -6500e3, 6500e3, 521
n_angles  = 100
n_pts     = 100
zposition = 2.5e3
energy_min = 1.0
energy_max = 100.0
energy_max = 40.0
osc_probs_nu, osc_probs_antinu = Neutrino_Oscillations.produce_neutrino_oscillation_probabilities(minX, maxX, nX, n_angles, n_pts, zposition, energy_min=energy_min, energy_max=energy_max)
# neutrinos, NMO
Pνe2νe = osc_probs_nu[:, :, 1, 1]' 
Pνe2νμ = osc_probs_nu[:, :, 1, 2]' 
Pνe2ντ = osc_probs_nu[:, :, 1, 3]' 
Pνμ2νe = osc_probs_nu[:, :, 2, 1]' 
Pνμ2νμ = osc_probs_nu[:, :, 2, 2]' 
Pνμ2ντ = osc_probs_nu[:, :, 2, 3]' 
# antineutrinos, NMO
Paνe2aνe = osc_probs_antinu[:, :, 1, 1]' 
Paνe2aνμ = osc_probs_antinu[:, :, 1, 2]' 
Paνe2aντ = osc_probs_antinu[:, :, 1, 3]' 
Paνμ2aνe = osc_probs_antinu[:, :, 2, 1]' 
Paνμ2aνμ = osc_probs_antinu[:, :, 2, 2]' 
Paνμ2aντ = osc_probs_antinu[:, :, 2, 3]' 

# plots for neutrinos, NMO
angles = range(-1, 0, length=100)
titles = (
    L"P($\nu_e \rightarrow \nu_e$)", L"P($\nu_e \rightarrow \nu_\mu$)", L"P($\nu_e \rightarrow \nu_\tau$)",
    L"P($\nu_\mu \rightarrow \nu_e$)", L"P($\nu_\mu \rightarrow \nu_\mu$)", L"P($\nu_\mu \rightarrow \nu_\tau$)"
)
all_probs = (
    Pνe2νe, Pνe2νμ, Pνe2ντ,
    Pνμ2νe, Pνμ2νμ, Pνμ2ντ
)
plot_list = []
for idx in 1:6
    # determine grid positions for labels
    is_bottom_row = idx > 3
    is_left_column = (idx == 1 || idx == 4)
    is_right_column = (idx ==3 || idx == 6)
    p = heatmap(
        energies, angles, all_probs[idx],              
        title = titles[idx],
        xlabel = is_bottom_row ? "Energy/GeV" : "",
        ylabel = is_left_column ? L"\cos\theta" : "",
        colorbar_title = is_right_column ? "Oscillation probability" : "",
        xscale = :log10,
        clim = (0, 1),
        cmap = :viridis
    )
    push!(plot_list, p)
end
final_plot = plot(
    plot_list..., 
    layout = (2, 3), 
    size = (1000, 400), 
    margin = 1Plots.mm,
)

# plots for antineutrinos, NMO
titles = (
    L"P($\bar{\nu}_e \rightarrow \bar{\nu}_e$)", L"P($\bar{\nu}_e \rightarrow \bar{\nu}_\mu$)", L"P($\bar{\nu}_e \rightarrow \bar{\nu}_\tau$)",
    L"P($\bar{\nu}_\mu \rightarrow \bar{\nu}_e$)", L"P($\bar{\nu}_\mu \rightarrow \bar{\nu}_\mu$)", L"P($\bar{\nu}_\mu \rightarrow \bar{nu}_\tau$)"
)
all_probs = (
    Paνe2aνe, Paνe2aνμ, Paνe2aντ,
    Paνμ2aνe, Paνμ2aνμ, Paνμ2aντ
)
plot_list = []
for idx in 1:6
    # determine grid positions for labels
    is_bottom_row = idx > 3
    is_left_column = (idx == 1 || idx == 4)
    is_right_column = (idx ==3 || idx == 6)
    p = heatmap(
        energies, angles, all_probs[idx],              
        title = titles[idx],
        xlabel = is_bottom_row ? "Energy/GeV" : "",
        ylabel = is_left_column ? L"\cos\theta" : "",
        colorbar_title = is_right_column ? "Oscillation probability" : "",
        xscale = :log10,
        clim = (0, 1),
        cmap = :viridis
    )
    push!(plot_list, p)
end
final_plot = plot(
    plot_list..., 
    layout = (2, 3), 
    size = (1000, 400), 
    margin = 1Plots.mm
)





energies = logrange(0.1, 100, 100)
energies = logrange(1, 40, 100)
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





# all together → to get interacting events histograms
energies = logrange(1, 40, 100)
ΔE = diff(energies)
# temporary fix
push!(ΔE, ΔE[end])
#pushfirst!(ΔE, ΔE[1])
cosθs = range(-1, 0, length=100)
Δcosθ = diff(cosθs)[1] 
# data-taking time in seconds, we use a year here
time  = 365.25 * 24.0 * 60.0 * 60.0 
# proton mass in kg
p_m   = 1.67262e-27
# detector mass in kg
det_m = 1.0e9 
# detector exposure in seconds * kg 
ξ_Det = time * det_m
# Interacting electron neutrino events
constant = (ξ_Det ./ p_m) .* 1.0e-42        .* 2.0 .* π .* Δcosθ      .*10
Edeppart = cs_νe_CC .* energies             .* (energies.^2) .*ΔE
int_evts_νe = constant .* Edeppart' .* (flux_νe .* Pνe2νe .+ flux_νμ .* Pνμ2νe)
Edeppartanti = cs_antiνe_CC .* energies             .* (energies.^2) .*ΔE
int_evts_antiνe = constant .* Edeppart' .* (flux_antiνe .* Paνe2aνe .+ flux_antiνμ .* Paνμ2aνe)
heatmap(energies, cosθs, int_evts_νe+int_evts_antiνe,
        xscale = :log10,
        xlabel = "Energy/GeV",
        ylabel = L"\cos\theta")
# Interacting tau neutrino events
Edeppart = cs_ντ_CC .* energies             .* (energies.^2) .*ΔE
int_evts_ντ = constant .* Edeppart' .* (flux_νe .* Pνe2ντ .+ flux_νμ .* Pνμ2ντ)
Edeppartanti = cs_antiντ_CC .* energies             .* (energies.^2) .*ΔE
int_evts_antiντ = constant .* Edeppart' .* (flux_antiνe .* Paνe2aντ .+ flux_antiνμ .* Paνμ2aντ)
heatmap(energies, cosθs, int_evts_ντ+int_evts_antiντ,
        xscale = :log10,
        xlabel = "Energy/GeV",
        ylabel = L"\cos\theta")












#matrix = [
#    1.0  2.0 ;  # Row 1
#    3.0  4.0    # Row 2
#]
#vector = [1.0, 2.0]
#vector .* matrix
#vector' .* matrix

