using CSV
using DataFrames

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
    raw_df = CSV.read(csv_path, DataFrame, header=false)

   # ensure the number of rows matches the expected grid size
    expected_rows = nEbins * nθbins
    if nrow(raw_df) != expected_rows
        error("Dimension mismatch: CSV has $(nrow(raw_df)) rows, but nEbins * nθbins = $expected_rows.")
    end

    # convert DataFrame to a 2D Matrix, then reshape to a 3D grid
    flux_3d = reshape(Matrix(raw_df), nEbins, nθbins, 5)
 
    # Slice the 3D grid to extract the 2D planes for each specific flavor
    nuflux_νe     = flux_3d[:, :, 4]
    nuflux_νμ     = flux_3d[:, :, 2]
    nuflux_antiνe = flux_3d[:, :, 5]
    nuflux_antiνμ = flux_3d[:, :, 3]

    return nuflux_νe, nuflux_νμ, nuflux_antiνe, nuflux_antiνμ
end



