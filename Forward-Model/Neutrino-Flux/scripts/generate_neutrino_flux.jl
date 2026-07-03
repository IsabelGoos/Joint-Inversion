module Neutrino_Flux
export read_neutrino_flux_table

using CSV
using DataFrames
using Interpolations

"""
    read_neutrino_flux_table(filename::String, nEbins::Int, nθbins::Int, has_header::Bool)

Read a neutrino flux CSV file from the `../data` directory 
and create a neutrino flux array for each neutrino flavor and type.
The amount of tau neutrinos and antineutrinos is negligible.

# Arguments
- `filename::String`: The name of the data file (without the `.csv` extension).
- `nEbins::Int`: The number of energy bins.
- `nθbins::Int`: The number of zenith angle bins.
- `has_header::Bool`: Is true if the dataset has a header. 

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
    intepolate_flux_at_bin_centers(energies::Vector{Float64}, flux::Matrix{Float64}, energies_in_log::Bool)

Calculates energy bin centers, interpolates the log(flux) at log(energies).

# Returns
A tuple containing
1. bin centers
2. interpolated flux values evaluated at those centers
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
    # evaluate the interpolator at the log of the bin centers, then exponentiate back
    interpolated_flux = exp.(hcat([interp.(log.(bin_centers)) for interp in interpolators]...))   
    return bin_centers, interpolated_flux

end

function produce_neutrino_flux_table(model::String)

    #if model = "Daemonflux" etc

    return Nothing

end

end