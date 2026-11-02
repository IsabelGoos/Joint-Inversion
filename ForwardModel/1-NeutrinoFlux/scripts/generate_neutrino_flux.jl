include("config.jl")

using Downloads
using CodecZlib
using DelimitedFiles
using Interpolations: linear_interpolation, Flat
using PythonCall

# ============================================================================
# Main function
# ============================================================================

"""
    produce_neutrino_flux(bin_centers_arrays::Tuple)
    produce_neutrino_flux(user_params::AbstractDict{String})

Produce the atmospheric neutrino flux on a given binning, using either the
default neutrino flux model and parameters or substituting with `user_params`.

# Arguments
- `bin_centers_arrays::Tuple`: binning `(energies, cosθs)`. The flux is
  computed with the current default model and parameters.
- `user_params::AbstractDict{String}`: dictionary with the following keys:
  - `bin_centers_arrays` (required): binning `(energies, cosθs)`.
  - `model` (optional): replaces the default model.
  - any other key (optional): a parameter of the chosen model that
    replaces its default value.

# Returns 
a 4-tuple (or 8-tuple, if `return_uncertainties::Bool` is true) with the following matrices in the shape (nθ, nE):
- `NuMu_flux` 
- `Nue_flux`
- `AntiNuMu_flux` 
- `AntiNue_flux`
(followed by the corresponding uncertainties if `return_uncertainties::Bool` is true, in the same order).
"""
produce_neutrino_flux(bin_centers_arrays::NTuple{2, AbstractVector{<:Real}}) =
    produce_neutrino_flux(Dict{String, Any}("bin_centers_arrays" => bin_centers_arrays))

function produce_neutrino_flux(user_params::Dict{String, Any})

    # Make sure that user_params contains the bin centers.
    haskey(user_params, "bin_centers_arrays") ||
        throw(ArgumentError("The key \"bin_centers_arrays\" is required: (energies, cosθs)."))

    # Pick the model and parameter set.
    if haskey(user_params, "model")
        model  = _validate_neutrino_flux_model(user_params["model"])
        params = MODEL_2_PARAMSTYPE[model]()      
    else
        model  = DEFAULT_NUFLUX_MODEL[]
        params = DEFAULT_NUFLUX_PARAMS[]          
    end

    # Update parameters. Any problem will be raised by _update_parameters.
    updates = Dict{Symbol, Any}(
        Symbol(k) => v for (k, v) in user_params if k ∉ ("model", "bin_centers_arrays")
    )  
    params_validated = _validate_parameters(_update_parameters(params; updates...))

    # Dictionary expected by the _produce_*_neutrino_flux functions.
    nuflux_params = Dict{String, Any}(
        String(f) => getfield(params_validated, f) for f in fieldnames(typeof(params_validated))
    )
    nuflux_params["model"]              = model
    nuflux_params["bin_centers_arrays"] = user_params["bin_centers_arrays"]

    for (key, value) in nuflux_params
        println(key, " => ", value)
    end

    # Produce and return atmospheric neutrino flux.
    if model === :Honda
        nuflux = _produce_honda_neutrino_flux(nuflux_params)
    elseif model === :Daemonflux
        nuflux = _produce_daemonflux_neutrino_flux(nuflux_params)
    else  # :Cheng
        nuflux = _produce_cheng_neutrino_flux(nuflux_params)
    end

    # Interpolate if the models are Honda or Cheng.
    if model in (:Honda, :Cheng)
        E_orig    = nuflux["Energies"]
        cosθ_orig = nuflux["cosθs"]
        bins_new  = user_params["bin_centers_arrays"]
        nuflux_new = Dict{String, Any}( 
            key => _get_interpolated_flux(nuflux[key], E_orig, cosθ_orig, model, bins_new)
            for key in FLUX_KEYS
        )
        nuflux_new["Energies"] = collect(bins_new[1])
        nuflux_new["cosθs"]    = collect(bins_new[2])

        nuflux = nuflux_new
    end

    return nuflux
end

# ============================================================================
# Individual functions for each neutrino flux model
# ============================================================================

"""
    _produce_honda_neutrino_flux(nuflux_params::Dict{String, Any})

Produce the atmospheric neutrino flux for all relevant neutrino types using the Honda model.
"""
function _produce_honda_neutrino_flux(nuflux_params::Dict{String, Any})

    # Read parameters and build url from where to fetch the data.
    url_honda  = _build_flux_url(nuflux_params)
    angles_sep = nuflux_params["angles_separation"]
    
    # Download.
    decompressed_stream = download_nuflux(url_honda)
    
    # Read the text grid into a matrix, omitting any lines with text.
    lines = readlines(decompressed_stream)
    close(decompressed_stream)
    numeric_lines = filter(line -> occursin(r"^\s*[0-9.-]", line), lines)
    clean_text_block = join(numeric_lines, "\n")
    matrix_2d = readdlm(IOBuffer(clean_text_block), Float64)
    
    # Reshape, based on angular separation mode.

    # :averaged_ϕ-grid layout: 5 columns x 101 energies x 20 cosθ steps
    if angles_sep === :averaged_ϕ

        n_E    = 101
        n_cosθ = 20
        flux_reshaped = reshape(matrix_2d, n_E, n_cosθ, 5)

        # Extract matrices.
        results = Dict{String, Any}(
            "NuMu_flux"     => collect(flux_reshaped[:, :, 2]),
            "Nue_flux"      => collect(flux_reshaped[:, :, 4]),
            "AntiNuMu_flux" => collect(flux_reshaped[:, :, 3]),
            "AntiNue_flux"  => collect(flux_reshaped[:, :, 5]),
            "Energies"      => flux_reshaped[:, 1, 1],
            "cosθs"         => range(-1.0, 1.0, length=n_cosθ)
        )

    elseif angles_sep === :variable_ϕ
        # TODO: Implement the reshaping logic for the azimuth-dependent grid
        error("angles_separation_hf mode :variable_ϕ is not yet implemented.")

    elseif angles_sep === :averaged_ϕθ
        # TODO: Implement the reshaping logic for the fully integrated/averaged grid
        error("angles_separation_hf mode :averaged_ϕθ is not yet implemented.")

    else
        error("Invalid angles_separation_hf option: $angles_sep. Use :variable_ϕ, :averaged_ϕ, or :averaged_ϕθ.")
    end

    return results
end

"""
    _produce_daemonflux_neutrino_flux(nuflux_params::Dict{String, Any})

Produce the atmospheric neutrino flux for all relevant neutrino types using the Daemonflux model.
"""
function _produce_daemonflux_neutrino_flux(nuflux_params::Dict{String, Any})

    # Read some inputs.
    E_new, cosθ_new      = nuflux_params["bin_centers_arrays"]
    flux_mode            = nuflux_params["flux_mode"]
    return_uncertainties = nuflux_params["return_uncertainties"]

    # Import and instantiate Python daemonflux class.
    daemonflux = pyimport("daemonflux")
    nuflux = daemonflux.Flux(location="IceCube", 
                             use_calibration=true, 
                             debug=1)

    # Convert to degrees, as needed by daemonflux.
    θ_new = string.(rad2deg.(acos.(cosθ_new)))

    # Pre-allocate output arrays.
    nE = length(E_new)
    nθ = length(θ_new)
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

    # Define key names depending on the chosen flux_mode. 
    prefix = flux_mode === :total ? "total_" : ""
    keys = (
        numu    = prefix * "numu",
        nue     = prefix * "nue",
        anumu   = prefix * "antinumu",
        anue    = prefix * "antinue" 
    )

    # Set daemonflux parameters.
    raw_params_df = get(nuflux_params, "params_df", nothing)
    params = (isnothing(raw_params_df) || isempty(raw_params_df)) ? 
             _set_dflux_params(nuflux) : 
             _set_dflux_params(nuflux, raw_params_df)

    # Energy scaling.
    E_scaling = 1.0e4 ./ (E_new .^ 3)

    # Fetch raw daemonflux values.
    @views for (i, θ) in enumerate(θ_new)
        # flux values
        NuMu_flux[i, :]     .= pyconvert(Vector{Float64}, nuflux.flux(E_new, θ, keys.numu,  params=params)) .* E_scaling
        Nue_flux[i, :]      .= pyconvert(Vector{Float64}, nuflux.flux(E_new, θ, keys.nue,   params=params)) .* E_scaling
        AntiNuMu_flux[i, :] .= pyconvert(Vector{Float64}, nuflux.flux(E_new, θ, keys.anumu, params=params)) .* E_scaling
        AntiNue_flux[i, :]  .= pyconvert(Vector{Float64}, nuflux.flux(E_new, θ, keys.anue,  params=params)) .* E_scaling
        # flux uncertainties
        # E_scaling relevant here too? TODO: check
        if return_uncertainties
            NuMu_flux_err[i, :]     .= pyconvert(Vector{Float64}, nuflux.error(E_new, θ, keys.numu))
            Nue_flux_err[i, :]      .= pyconvert(Vector{Float64}, nuflux.error(E_new, θ, keys.nue))
            AntiNuMu_flux_err[i, :] .= pyconvert(Vector{Float64}, nuflux.error(E_new, θ, keys.anumu))
            AntiNue_flux_err[i, :]  .= pyconvert(Vector{Float64}, nuflux.error(E_new, θ, keys.anue))
        end
    end

    # Prepare dictionary for output.
    results = Dict{String, Any}(
        "NuMu_flux"     => NuMu_flux,
        "Nue_flux"      => Nue_flux,
        "AntiNuMu_flux" => AntiNuMu_flux,
        "AntiNue_flux"  => AntiNue_flux
    )
    if return_uncertainties
        results["NuMu_flux_err"]     = NuMu_flux_err
        results["Nue_flux_err"]      = Nue_flux_err
        results["AntiNuMu_flux_err"] = AntiNuMu_flux_err
        results["AntiNue_flux_err"]  = AntiNue_flux_err
    end
    return results

end

"""
    _set_dflux_params(flux_obj)

Initialize a dictionary of Daemonflux parameters set to their default values.

# Arguments
- `flux_obj`: A Python object representing a Daemonflux instance.
"""
function _set_dflux_params(flux_obj)

    # Fetch all the Daemonflux parameter names. 
    params = pyconvert(Vector{String}, flux_obj.params.known_parameters)

    # Set all parameters to 0 
    # (meaning that their difference with respect to the default parameters is zero).
    params_dict = Dict(p => 0.0 for p in params) 
    
    return params_dict

end

"""
    _set_dflux_params(flux_obj, params_dict)

Initialize a dictionary of Daemonflux parameters to those values given in `params_dict`.
Those parameters not mentioned there are set to default values.

# Arguments
- `flux_obj`: A Python object representing a Daemonflux instance.
- `params_dict`: A dictionary with those pairs of parameter names and values which should be modified.
The possible parameter names are K+_158G, K+_2P, K+_31G, K-_158G, K-_2P, K-_31G, n_158G, n_2P, p_158G, p_2P, 
pi+_158G, pi+_20T, pi+_2P, pi+_31G, pi-_158G, pi-_20T, pi-_2P, pi-_31G, GSF_1, GSF_2, GSF_3, GSF_4, GSF_5, GSF_6.
The get_dflux_param_value function helps find specific values.

# Returns
- `Dict{String, Float64}` 
"""
function _set_dflux_params(flux_obj, params_df_dict::Union{Nothing, Dict{String, Any}})

    # Set default Daemonflux parameters.
    full_params = _set_dflux_params(flux_obj)

    # If params_df_dict is nothing, just return the defaults without trying to merge.
    if isnothing(params_df_dict)
        return full_params
    end

    # Otherwise, update with pairs from params_dict.
    merge!(full_params, params_df_dict)

    return full_params

end

"""
    _produce_cheng_neutrino_flux(nuflux_params::Dict{String, Any})

Produce the atmospheric neutrino flux for all relevant neutrino types using the Cheng model.
"""
function _produce_cheng_neutrino_flux(nuflux_params::Dict{String, Any})

    # Read parameters and build url from where to fetch the data.
    url_cheng  = _build_flux_url(nuflux_params)
    angles_sep = nuflux_params["angles_separation"]

    # Download.
    decompressed_stream = download_nuflux(url_cheng)

    # Read the text grid into a matrix, omitting any lines with text.
    lines = readlines(decompressed_stream)
    close(decompressed_stream)
    numeric_lines = filter(line -> occursin(r"^\s*[0-9.-]", line), lines)
    clean_text_block = join(numeric_lines, "\n")
    matrix_2d = readdlm(IOBuffer(clean_text_block), Float64)

    # Reshape, based on angular separation mode.

    # :averaged_ϕ-grid layout: 5 columns x 101 energies x 20 cosθ steps
    if angles_sep === :averaged_ϕ

        n_E    = 121
        n_cosθ = 20
        flux_reshaped = reshape(matrix_2d, n_E, n_cosθ, 5)

        # Extract matrices.
        results = Dict{String, Any}(
            "NuMu_flux"     => collect(flux_reshaped[:, :, 2]),
            "Nue_flux"      => collect(flux_reshaped[:, :, 4]),
            "AntiNuMu_flux" => collect(flux_reshaped[:, :, 3]),
            "AntiNue_flux"  => collect(flux_reshaped[:, :, 5]),
            "Energies"      => flux_reshaped[:, 1, 1],
            "cosθs"         => range(-1.0, 1.0, length=n_cosθ)
        )

    elseif angles_sep === :variable_ϕ
        # TODO: Implement the reshaping logic for the azimuth-dependent grid
        error("angles_separation_hf mode :variable_ϕ is not yet implemented.")

    elseif angles_sep === :averaged_ϕθ
        # TODO: Implement the reshaping logic for the fully integrated/averaged grid
        error("angles_separation_hf mode :averaged_ϕθ is not yet implemented.")

    else
        error("Invalid angles_separation_hf option: $angles_sep. Use :variable_ϕ, :averaged_ϕ, or :averaged_ϕθ.")
    end

    return results
end

# ============================================================================
# Function that builds the urls of the different neutrino flux models
# ============================================================================

"""
    _params_from_dict(nuflux_params::Dict{String, Any})

Rebuild the parameter struct from the dictionary `nuflux_params`.
"""
function _params_from_dict(nuflux_params::Dict{String, Any})
    ParamsType = MODEL_2_PARAMSTYPE[nuflux_params["model"]]
    fields = (Symbol(k) => v for (k, v) in nuflux_params
              if k ∉ ("model", "bin_centers_arrays"))
    return ParamsType(; fields...)
end

_build_flux_url(nuflux_params::Dict{String, Any}) =
    _build_flux_url(_params_from_dict(nuflux_params))

function _build_flux_url(p::HondaParameters)

    loc   = LOCATION_MAPPING_H[p.location]
    seas  = SEASON_MAPPING[p.season]
    ang   = ANGLES_MAPPING[p.angles_separation]
    sol   = SOLAR_MAPPING[p.solar_activity].honda 
    mtn   = MOUNTAIN_MAPPING[p.mountain_overburden]

    return "$HONDA_BASE_URL/$loc-$seas-$ang$mtn-$sol.d.gz"
end

function _build_flux_url(p::ChengParameters)

    loc1  = LOCATION_MAPPING_C[p.location][1]
    loc2  = LOCATION_MAPPING_C[p.location][2]
    seas  = SEASON_MAPPING[p.season]
    ang   = ANGLES_MAPPING[p.angles_separation]
    sol1  = SOLAR_MAPPING[p.solar_activity].cheng 
    sol2  = SOLAR_MAPPING[p.solar_activity].honda 
    mu    = MUON_PROPA_MAPPING[p.muon_propa_earth]

    return "$CHENG_BASE_URL/$loc1/$(sol1)_$mu/$loc2-$seas-$ang-$sol2.d"
end

"""
    download_nuflux(url::String)

Download a neutrino flux and return a stream of its decompressed contents.
"""
function download_nuflux(url::String)
    
    # Download into an in-memory buffer and get the raw bytes.
    io_buffer = IOBuffer()
    try
        Downloads.download(url, io_buffer; timeout = 60)
    catch e 
        if e isa Downloads.RequestError
            error("Failed to download the flux table. The URL does not exist or the server is down.\n" *
                  "URL attempted: $url\nDetails: $(sprint(showerror, e))")        
        else
            rethrow(e) 
        end
    end
    data = take!(io_buffer)
    
    # Gzip files start with the bytes 0x1f 0x8b.
    is_gzip = length(data) >= 2 && data[1] == 0x1f && data[2] == 0x8b

    # Decompress, if possible.
    decompressed_stream = is_gzip ? GzipDecompressorStream(IOBuffer(data)) : IOBuffer(data)

    return decompressed_stream
end

# ============================================================================
# Functions to interpolate neutrino fluxes using the provided bins
# ============================================================================

"""
    _get_interpolated_flux(Phi_orig, E_orig, cosθ_orig, model, bin_centers_arrays)

Interpolate a flux table onto a new binning.

`Phi_orig` is a matrix defined on the original grid `(E_orig, cosθ_orig)`.
It is interpolated onto the grid `(E_new, cosθ_new)` given by `bin_centers_arrays`.

The interpolation is linear in `log(E)` and `log(Phi)` and linear in `cosθ`.

# Arguments
- `Phi_orig::AbstractMatrix{<:Real}`: flux table.
- `E_orig::AbstractVector{<:Real}`: energies of the original table, in GeV.
- `cosθ_orig::AbstractVector{<:Real}`: cosθ values of the original table
  (ascending or descending, in the same order as the columns of `Phi_orig`).
- `model::Symbol`: atmospheric neutrino flux model used.
- `bin_centers_arrays::NTuple{2, AbstractVector{<:Real}}`: new binning `(E_new, cosθ_new)`.

# Returns
A matrix of size `length(E_new) × length(cosθ_new)`.
"""
function _get_interpolated_flux(Phi_orig::AbstractMatrix{<:Real}, 
                                E_orig::AbstractVector{<:Real}, 
                                cosθ_orig::AbstractVector{<:Real}, 
                                model::Symbol, 
                                bin_centers_arrays::NTuple{2, AbstractVector{<:Real}})

    # Check if the new bin centers are within the allowed range.
    E_range, cosθ_range = MODEL_2_RANGETYPE[model]
    E_new,   cosθ_new   = bin_centers_arrays
    _check_bins_in_range(E_new, cosθ_new, E_range, cosθ_range)

    # The interpolation nodes must be in ascending order.
    order     = sortperm(cosθ_orig)
    cosθ_orig = cosθ_orig[order]
    Phi_orig  = Phi_orig[:, order]

    # Interpolate in log-log space for energy and flux. Atmospheric neutrino
    # spectra fall roughly as a power law in terms of energy. 
    logE_orig   = log.(E_orig)
    logPhi_orig = log.(Phi_orig)
    itp = linear_interpolation((logE_orig, cosθ_orig), logPhi_orig;
                               extrapolation_bc = Flat())

    # Evaluate on the new grid and go back from log(Phi) to Phi.
    Phi_new = [exp(itp(log(E), cosθ)) for E in E_new, cosθ in cosθ_new]
    return Phi_new
end

"""
    _check_bins_in_range(E_new, cosθ_new, E_range, cosθ_range; rtol=1e-8)

Check that a requested binning lies inside the allowed range, 
so that interpolation never needs to extrapolate.

# Arguments
- `E_new`: requested energies, in GeV.
- `cosθ_new`: requested cosines of the zenith angle.
- `E_range`: allowed energy range.
- `cosθ_range`: allowed angular range.

# Keywords
- `rtol::Real = 1e-8`: tolerance.
"""
function _check_bins_in_range(E_new, cosθ_new, E_range, cosθ_range; rtol=1e-8)

    # Accepted energy and angular ranges.
    E_lo,  E_hi  = E_range
    ct_lo, ct_hi = cosθ_range

    # Proposed energy and angular ranges.
    Emin,  Emax  = extrema(E_new)
    ctmin, ctmax = extrema(cosθ_new)

    # Check, with some tolerance, if the proposed ranges are within the accepted ones.
    (Emin >= E_lo * (1 - rtol) && Emax <= E_hi * (1 + rtol)) ||
        throw(ArgumentError("Energies [$Emin, $Emax] GeV are outside the allowed range [$E_lo, $E_hi] GeV."))
    (ctmin >= ct_lo - rtol && ctmax <= ct_hi + rtol) ||
        throw(ArgumentError("cosθ values [$ctmin, $ctmax] are outside the allowed range [$ct_lo, $ct_hi]."))
    return nothing
end


