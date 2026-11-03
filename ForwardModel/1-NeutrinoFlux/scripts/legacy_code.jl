"""
    read_neutrino_flux_table(path; nEbins, nθbins, has_header=false,
        cosθgrid=range(-1, 1, length=nθbins), energies_in_log=true)

Read a legacy five-column CSV table and return a compact `NamedTuple` with
`native` and energy-bin-centred `interpolated` tables. Flux matrices always use
the flexOPT convention `(energy, cosθ)`.

`path` may be an explicit CSV path, a filename relative to
`dataInput/neutrinoFluxData`, or a basename without `.csv`.
"""
function read_neutrino_flux_table(
    path::AbstractString;
    nEbins::Integer,
    nθbins::Integer,
    has_header::Bool=false,
    cosθgrid=range(-1.0, 1.0, length=nθbins),
    energies_in_log::Bool=true,
)
    csv_path = _resolve_neutrino_flux_csv(path)
    raw_df = CSV.read(csv_path, DataFrame; header=has_header)
    expected_rows = nEbins * nθbins
    nrow(raw_df) == expected_rows ||
        throw(DimensionMismatch(
            "CSV has $(nrow(raw_df)) rows; expected $nEbins × $nθbins = $expected_rows",
        ))
    ncol(raw_df) == 5 ||
        throw(DimensionMismatch(
            "legacy neutrino-flux CSV must contain five columns; got $(ncol(raw_df))",
        ))

    flux_3d = reshape(Float64.(Matrix(raw_df)), nEbins, nθbins, 5)
    energies = collect(flux_3d[:, 1, 1])
    cosines = Float64.(collect(cosθgrid))
    length(cosines) == nθbins ||
        throw(DimensionMismatch("cosθgrid must contain $nθbins values"))

    native = (
        energies=energies,
        cosθgrid=cosines,
        νe=collect(flux_3d[:, :, 4]),
        νμ=collect(flux_3d[:, :, 2]),
        antiνe=collect(flux_3d[:, :, 5]),
        antiνμ=collect(flux_3d[:, :, 3]),
    )
    interpolated_fluxes = map((native.νe, native.νμ, native.antiνe, native.antiνμ)) do flux
        _, values = interpolate_flux_at_bin_centers(
            energies,
            flux,
            energies_in_log,
        )
        values
    end
    bin_centers, _ = interpolate_flux_at_bin_centers(
        energies,
        native.νe,
        energies_in_log,
    )
    interpolated = (
        energies=bin_centers,
        cosθgrid=cosines,
        νe=interpolated_fluxes[1],
        νμ=interpolated_fluxes[2],
        antiνe=interpolated_fluxes[3],
        antiνμ=interpolated_fluxes[4],
    )
    return (; native, interpolated)
end

function _resolve_neutrino_flux_csv(path::AbstractString)
    candidates = String[path]
    endswith(lowercase(path), ".csv") || push!(candidates, path * ".csv")
    data_dir = normpath(joinpath(@__DIR__, "..", "..", "dataInput", "neutrinoFluxData"))
    append!(candidates, joinpath.(Ref(data_dir), basename.(copy(candidates))))
    match = findfirst(isfile, candidates)
    isnothing(match) &&
        throw(ArgumentError(
            "neutrino-flux CSV not found; tried $(join(candidates, ", "))",
        ))
    return abspath(candidates[match])
end

# Positional compatibility without the former ten-element return value.
function read_neutrino_flux_table(
    path::AbstractString,
    nEbins::Integer,
    nθbins::Integer,
    has_header::Bool,
)
    return read_neutrino_flux_table(
        path;
        nEbins=nEbins,
        nθbins=nθbins,
        has_header=has_header,
    )
end
