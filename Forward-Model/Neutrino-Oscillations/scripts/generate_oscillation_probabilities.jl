module Neutrino_Oscillations
export set_oscillation_parameters, produce_neutrino_oscillation_probabilities

const FLEXOPT_DIR = "/Users/igoos/Desktop/projects/flexOPT" # TO-DO: adapt this path to FlexOPT
include(joinpath(FLEXOPT_DIR, "src", "Neurthino.jl"))
using .Neurthino
include(joinpath(FLEXOPT_DIR, "src", "NeurthinoBack", "Matter.jl"))
include(joinpath(FLEXOPT_DIR, "src", "NeurthinoBack", "premFunctions.jl"))
using Interpolations
ParamFile = "../config/testparam.csv" # TO-DO: adapt this path, in case it changes 
include("generate_earth_model.jl")


"""
    set_oscillation_parameters(ordering=:normal)

Initializes a 3-flavor neutrino oscillation object using the `Neurthino` framework
and the latest global best-fit values from NuFIT v6.1.
Accepts either `:normal` or `:inverted` for the `ordering` keyword.

### Examples:
set_oscillation_parameters()
set_oscillation_parameters(:inverted)
"""
function set_oscillation_parameters(ordering::Symbol=:normal)
 
    # initialize a 3-flavor oscillation object
    osc = Neurthino.OscillationParameters(3)
    # NuFIT v6.1 normal ordering with SK atmospheric data
    if ordering === :normal
        # mixing angles (radians)
        Neurthino.setθ!(osc, 1=>2, 0.589)   # θ12 ≈ 33.76°
        Neurthino.setθ!(osc, 1=>3, 0.150)   # θ13 ≈ 8.62°
        Neurthino.setθ!(osc, 2=>3, 0.756)   # θ23 ≈ 43.29° 
        # CP-violating phase (radians)
        Neurthino.setδ!(osc, 1=>3, 3.700)   # δCP ≈ 212°
        # Mass-squared splittings (eV²)
        Neurthino.setΔm²!(osc, 2=>1, 7.537e-5)   
        Neurthino.setΔm²!(osc, 3=>2, 2.511e-3)  
        
    elseif ordering === :inverted
        # Mixing angles (radians)
        Neurthino.setθ!(osc, 1=>2, 0.589)   # θ12 ≈ 33.76°
        Neurthino.setθ!(osc, 1=>3, 0.151)   # θ13 ≈ 8.65°
        Neurthino.setθ!(osc, 2=>3, 0.836)   # θ23 ≈ 47.90° 
        # CP-violating phase (radians)
        Neurthino.setδ!(osc, 1=>3, 4.782)   # δCP ≈ 274°
        # Mass-squared splittings (eV²)
        Neurthino.setΔm²!(osc, 2=>1,  7.537e-5)  
        Neurthino.setΔm²!(osc, 3=>2, -2.483e-3) 
        
    else
        error("Unknown mass ordering: :\$ordering. Use :normal or :inverted.")

    end
    
    return osc
end

"""
    set_oscillation_parameters(θ12, θ13, θ23, δCP, m12, m23)

Initializes a 3-flavor neutrino oscillation object using the `Neurthino` framework
and user-defined oscillation parameters.

### Note on units:
* Mixing angles (`θ12`, `θ13`, `θ23`) & CP-violating phase (`δCP`): must be provided in degrees. 
* Mass-squared splittings (`m12`, `m23`): must be provided in eV^2.

### Example:
set_oscillation_parameters(34.0, 9.0, 43.0, 212.0, 7.4e-5, 2.5e-3)
"""
function set_oscillation_parameters(θ12, θ13, θ23, δCP, m12, m23)
 
    # initialize a 3-flavor oscillation object
    osc = Neurthino.OscillationParameters(3)

    # mixing angles (radians)
    Neurthino.setθ!(osc, 1=>2, θ12*pi/180)   
    Neurthino.setθ!(osc, 1=>3, θ13*pi/180)   
    Neurthino.setθ!(osc, 2=>3, θ23*pi/180)  
    # CP-violating phase (radians)
    Neurthino.setδ!(osc, 1=>3, δCP)  
    # Mass-squared splittings (eV²)
    Neurthino.setΔm²!(osc, 2=>1, m12)   
    Neurthino.setΔm²!(osc, 3=>2, m23)  

    return osc
end

"""
    set_oscillation_parameters(; ordering=:normal, kwargs...)

Configure a 3-flavor neutrino oscillation object using default oscillation parameters
and doing some manual modifications.

### Arguments:
* `ordering::Symbol`: The baseline mass ordering. Options are `:normal` (default) or `:inverted`.

### Keyword arguments (optional modifications):
* `θ12`, `θ13`, `θ23`: in degrees.
* `δCP`: CP-violating phase in degrees.
* `m12`, `m23`: Mass-squared splittings in eV^2.

### Example:
# Use inverted ordering as the baseline, but modify just the CP-violating phase:
set_oscillation_parameters(ordering=:inverted, δCP=1.5)
"""
function set_oscillation_parameters(; ordering::Symbol=:normal,
    θ12 = nothing,
    θ13 = nothing,
    θ23 = nothing,
    δCP = nothing,
    m12 = nothing,
    m23 = nothing)

    # initialize a 3-flavor oscillation object 
    # and assign it the oscillation parameters of the chosen ordering
    osc = set_oscillation_parameters(ordering)

    # modify the parameters of interest
    if !isnothing(θ12) Neurthino.setθ!(osc, 1=>2, θ12*pi/180) end
    if !isnothing(θ13) Neurthino.setθ!(osc, 1=>3, θ13*pi/180) end 
    if !isnothing(θ23) Neurthino.setθ!(osc, 2=>3, θ23*pi/180) end 
    if !isnothing(δCP) Neurthino.setδ!(osc, 1=>3, δCP) end 
    if !isnothing(m12) Neurthino.setΔm²!(osc, 2=>1, m12) end
    if !isnothing(m23) Neurthino.setΔm²!(osc, 3=>2, m23) end 

    return osc
end

"""
    produce_neutrino_oscillation_probabilities(minX, maxX, nX, n_angles, n_pts, zposition; kwargs...)

Generate neutrino and anti-neutrino oscillation probabilities across a 2D interpolated 
Earth density model for a range of log-spaced energies.

### Arguments:
- `minX`, `maxX`: Spatial bounds of the Cartesian model grid (in meters).
- `nX`: Number of grid nodes along the X axis.
- `n_angles`: Number of path trajectories (vectors) radiating from the detector.
- `n_pts`: Number of sampling points along each baseline trajectory.
- `zposition`: Depth of the detector below the Earth's surface (in meters).

### Optional Keyword Arguments:
- `ordering`: Subtype framework (`:normal` or `:inverted`).
- `energy_min`, `energy_max`: Range boundaries for the neutrino energy spectrum (in GeV).
- `n_energies`: Total number of log-spaced energy bins.

### Returns:
- `osc_probs_nu`: 4D array of probabilities for neutrinos.
- `osc_probs_antinu`: 4D array of probabilities for anti-neutrinos.
"""
function produce_neutrino_oscillation_probabilities(
    minX, maxX, nX, n_angles, n_pts, zposition;
    ordering::Symbol = :normal,
    energy_min = 1.0,
    energy_max = 100.0,
    n_energies = 100)

    # load  Earth density model
    interpolated_density, _ = Earth_Model.load_first_interpolated_density_model(minX, maxX, nX)
    
    # calculate spatial resolution 
    dR = (maxX - minX) / (nX - 1)
    
    # create neutrino paths
    paths = PREMcreationPaths(interpolated_density, n_angles, n_pts, zposition, dR)

    # set oscillation parameters
    osc = set_oscillation_parameters(ordering=ordering)
    
    # generate log-spaced energy vector 
    log_min = log10(energy_min)
    log_max = log10(energy_max)
    energies = 10 .^ range(log_min, log_max, length=n_energies)

    # compute neutrino oscillations 
    osc_probs_nu     = Neurthino.Pνν(osc, energies, paths, anti=false)
    osc_probs_antinu = Neurthino.Pνν(osc, energies, paths, anti=true)
    
    return osc_probs_nu, osc_probs_antinu
end

end