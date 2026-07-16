module Neutrino_Flux
export read_neutrino_flux_table, produce_neutrino_flux

using CSV
using DataFrames
using Interpolations
using PythonCall

"""
    read_neutrino_flux_table(filename::String, nEbins::Int, nθbins::Int, has_header::Bool)

Read a neutrino flux CSV file from the `../data` directory 
and create a neutrino flux array for each neutrino flavor and type.
The amount of tau neutrinos and antineutrinos is negligible.

# Arguments
- `filename`: The name of the data file (without the `.csv` extension).
- `nEbins`: The number of energy bins.
- `nθbins`: The number of zenith angle bins.
- `has_header`: Is true if the dataset has a header. 

# Returns
A tuple of four `nEbins` × `nθbins` Matrices representing the flux for:
1. `νe` - Electron neutrino
2. `νμ` - Muon neutrino
3. `antiνe` - Electron antineutrino
4. `antiνμ` - Muon antineutrino
"""
function read_neutrino_flux_table(filename::String, nEbins::Int, nθbins::Int, has_header::Bool)

    # construct path
    data_dir = joinpath(@__DIR__, "..", "data")
    csv_path = joinpath(data_dir, "$(filename).csv")

    # read the data (has_header=false ensures the first row isn't taken to be the description)
    raw_df = CSV.read(csv_path, DataFrame, header=has_header)

   # ensure the number of rows matches the expected grid size
    expected_rows = nEbins * nθbins
    if nrow(raw_df) != expected_rows
        error("Dimension mismatch: CSV has $(nrow(raw_df)) rows, but nEbins * nθbins = $expected_rows.")
    end

    # convert DataFrame to a 2D Matrix, then reshape to a 3D grid
    flux_3d = reshape(Matrix(raw_df), nEbins, nθbins, 5)
 
    # Slice the 3D grid to extract the 2D planes for each specific flavor
    energies      = flux_3d[:, 1, 1]
    nuflux_νμ     = flux_3d[:, :, 2]
    nuflux_antiνμ = flux_3d[:, :, 3]
    nuflux_νe     = flux_3d[:, :, 4]
    nuflux_antiνe = flux_3d[:, :, 5]

    # Compute values at bin centers
    bin_centers, nuflux_νμ_interp = interpolate_flux_at_bin_centers(energies, nuflux_νμ, true)
    -, nuflux_antiνμ_interp = interpolate_flux_at_bin_centers(energies, nuflux_antiνμ, true)
    -, nuflux_νe_interp     = interpolate_flux_at_bin_centers(energies, nuflux_νe, true)
    -, nuflux_antiνe_interp = interpolate_flux_at_bin_centers(energies, nuflux_antiνe, true)

    return bin_centers, nuflux_νe_interp, nuflux_νμ_interp, nuflux_antiνe_interp, nuflux_antiνμ_interp, energies, nuflux_νe, nuflux_νμ, nuflux_antiνe, nuflux_antiνμ

end

"""
    interpolate_flux_at_bin_centers(energies::Vector{Float64}, flux::Matrix{Float64}, energies_in_log::Bool)

Calculates energy bin centers, interpolates the log(flux) at log(energies).

# Arguments
- `energies`: The bin edges for the energy axis.
- `flux`: The flux data matrix.
- `energies_in_log`: If `true`, it assumes the input `energies` are evenly spaced in log-scale 
  andit calculates geometric centers. Otherwise, it calculates arithmetic centers.

# Returns
A tuple containing:
1. `bin_centers`: The calculated centers of the energy bins.
2. `interpolated_flux`: The flux evaluated at the bin centers.
"""
function interpolate_flux_at_bin_centers(energies::Vector{Float64}, flux::Matrix{Float64}, energies_in_log::Bool)
    
    # create 1D coordinate array for the y axis
    angles_indices = 1.0:size(flux, 2)

    # calculate the bin centers    
    if energies_in_log == true
        bin_centers = sqrt.(energies[1:end-1] .* energies[2:end]) 
    else
        bin_centers = 0.5 .* (energies[1:end-1] .+ energies[2:end])
    end

    # we interpolate log(flux) vs log(energies)
    interpolators = [interpolate(log.(energies), log.(col), FritschCarlsonMonotonicInterpolation()) for col in eachcol(flux)]
    # evaluate the interpolators at the log of the bin centers, then exponentiate back
    interpolated_flux = exp.(hcat([interp.(log.(bin_centers)) for interp in interpolators]...))   
    return bin_centers, interpolated_flux

end

"""
    produce_neutrino_flux(bin_centers_arrays; model, flux_mode)

Produce the atmospheric neutrino flux for all relevant neutrino types following a given `model`.

If the `model` `Daemonflux` is chosen, compute the atmospheric neutrino fluxes at `IceCube` 
using the Daemonflux model (arXiv:2303.00022).
`Icecube` is, for the moment, the only location for which upgoing neutrinos are computed (July 2026).
`use_calibration` is set to true so that the Daemonflux calibration is used.
If the `model` `Honda` is chosen, fetch the atmospheric neutrino fluxes from 
<http://www-rccn.icrr.u-tokyo.ac.jp/mhonda/public/nflx2014/index.html> (arXiv:1502.03916).
Neither model provides tau-neutrino or -antineutrino fluxes since they are negligible.

# Arguments
- `bin_centers_arrays::Tuple{AbstractArray, AbstractArray}`: a 2-tuple containing the bin centers 
for energy and cosθ. The bin centers for cosθ have to be given in increasing order.
- `model`: `Daemonflux` (default) or `Honda`.
- `flux_mode`: `conventional` or `total` (default). `total` means conventional + prompt (only the 
conventional componennt is calibrated).
- `return_uncertainties::Bool`: if true, return neutrino flux uncertainties. Default is false. 
- `df_params::Union{Nothing, Dict{String, <:Number}}`: Daemonflux parameters. By default all Daemonflux parameters
are set to their default values. Only those provided in the `df_params` are updated to the desired values.

# Returns 
a 4-tuple (or 8-tuple, if `return_uncertainties::Bool` is true) with the following matrices 
in the shape (nθ, nE):
- `NuMu_flux` 
- `Nue_flux`
- `AntiNuMu_flux` 
- `AntiNue_flux`
(followed by the corresponding uncertanties if `return_uncertainties::Bool` is true, in the same order).
"""
function produce_neutrino_flux(
    bin_centers_arrays::Tuple{AbstractArray, AbstractArray}; 
    model::Symbol=:Daemonflux, 
    flux_mode::Symbol=:total,
    return_uncertainties::Bool=false
    df_params::Union{Nothing, Dict{String, <:Number}} = nothing
)

    if (model === :Daemonflux)

        # Import and instantiate Python Flux class
        daemonflux = pyimport("daemonflux")
        flux = daemonflux.Flux(location="IceCube", 
                               use_calibration=true, 
                               debug=1)

        # Read bin centers for energy and cosθ
        Ebin_centers, cosθbin_centers = bin_centers_arrays
        # Convert to degrees
        θbin_centers = string.(rad2deg.(acos.(cosθbin_centers)))

        # Pre-allocate output arrays
        nE = length(Ebin_centers)
        nθ = length(θbin_centers)
        NuMu_flux     = zeros(Float64, nθ, nE)
        Nue_flux      = zeros(Float64, nθ, nE)
        AntiNuMu_flux = zeros(Float64, nθ, nE)
        AntiNue_flux  = zeros(Float64, nθ, nE)
        if return_uncertainties
            NuMu_flux_err     = zeros(Float64, nθ, nE)
            Nue_flux_err      = zeros(Float64, nθ, nE)
            AntiNuMu_flux_err = zeros(Float64, nθ, nE)
            AntiNue_flux_err  = zeros(Float64, nθ, nE)
        end

        # Adapt the key names depending on the chosen flux_mode 
        if flux_mode !== :total && flux_mode !== :conventional
            error("flux_mode $flux_mode is not a valid option. Use :total or :conventional.")
        end
        prefix = flux_mode === :total ? "total_" : ""
        keys = (
            numu    = prefix * "numu",
            nue     = prefix * "nue",
            anumu   = prefix * "antinumu",
            anue    = prefix * "antinue" 
        )

        # set Daemonflux parameters
        params = set_daemonflux_parameters(df_params)

        # Energy scaling
        E_scaling = 1e4 ./ (Ebin_centers .^ 3)
        # Fetch raw Daemonflux values
        for (i, θ) in enumerate(θbin_centers)
            # flux values
            raw_numu  = pyconvert(Vector{Float64}, flux.flux(Ebin_centers, θ, keys.numu,  params=params))
            raw_nue   = pyconvert(Vector{Float64}, flux.flux(Ebin_centers, θ, keys.nue,   params=params))
            raw_anumu = pyconvert(Vector{Float64}, flux.flux(Ebin_centers, θ, keys.anumu, params=params))
            raw_anue  = pyconvert(Vector{Float64}, flux.flux(Ebin_centers, θ, keys.anue,  params=params))
            NuMu_flux[i, :]     .= raw_numu  .* E_scaling
            Nue_flux[i, :]      .= raw_nue   .* E_scaling
            AntiNuMu_flux[i, :] .= raw_anumu .* E_scaling
            AntiNue_flux[i, :]  .= raw_anue  .* E_scaling
            if return_uncertainties
                # flux uncertainties
                NuMu_flux_err[i, :]     .= pyconvert(Vector{Float64}, flux.error(Ebin_centers, θ, keys.numu))
                Nue_flux_err[i, :]      .= pyconvert(Vector{Float64}, flux.error(Ebin_centers, θ, keys.nue))
                AntiNuMu_flux_err[i, :] .= pyconvert(Vector{Float64}, flux.error(Ebin_centers, θ, keys.anumu))
                AntiNue_flux_err[i, :]  .= pyconvert(Vector{Float64}, flux.error(Ebin_centers, θ, keys.anue))
            end
        end

    elseif (model === :Honda)
        error("model $model is not YET supported - work in progress.")
    else 
        error("model $model is not supported. Use :Daemonflux or :Honda.")
    end

    if return_uncertainties
        return (NuMu_flux,     Nue_flux,     AntiNuMu_flux,     AntiNue_flux), 
               (NuMu_flux_err, Nue_flux_err, AntiNuMu_flux_err, AntiNue_flux_err)
    else
        return NuMu_flux, Nue_flux, AntiNuMu_flux, AntiNue_flux
    end

end

end