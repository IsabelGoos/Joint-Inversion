const FLEXOPT_DIR = "/Users/igoos/Desktop/projects/flexOPT" # TO-DO: adapt this path to FlexOPT
include(joinpath(FLEXOPT_DIR, "src", "Neurthino.jl"))

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
        Neurthino.setΔm²!(osc, 1=>2, 7.537e-5)   
        Neurthino.setΔm²!(osc, 2=>3, 2.511e-3)  
        
    elseif ordering === :inverted
        # Mixing angles (radians)
        Neurthino.setθ!(osc, 1=>2, 0.589)   # θ12 ≈ 33.76°
        Neurthino.setθ!(osc, 1=>3, 0.151)   # θ13 ≈ 8.65°
        Neurthino.setθ!(osc, 2=>3, 0.836)   # θ23 ≈ 47.90° 
        # CP-violating phase (radians)
        Neurthino.setδ!(osc, 1=>3, 4.782)   # δCP ≈ 274°
        # Mass-squared splittings (eV²)
        Neurthino.setΔm²!(osc, 1=>2,  7.537e-5)  
        Neurthino.setΔm²!(osc, 2=>3, -2.483e-3) 
        
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
    Neurthino.setΔm²!(osc, 1=>2, m12)   
    Neurthino.setΔm²!(osc, 2=>3, m23)  

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
function set_oscillation_parameters(; ordering::Symbol = :normal,
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
    if !isnothing(m12) Neurthino.setΔm²!(osc, 1=>2, m12) end
    if !isnothing(m23) Neurthino.setΔm²!(osc, 2=>3, m23) end 

    return osc
end






function linkWithNeurthinoPREM_YD(cos_θ,energies)
  
    paths = PREMcreationPaths(n_vectors, zposition)
    matprobs2 = Neurthino.Pνν(osc, energies, paths)
    matprobs2anti = Neurthino.Pνν(osc, energies, paths,anti=true)
    return matprobs2, matprobs2anti
end

"""
    generate_oscillation_probabilitiese(filename::String, nEbins::Int, nθbins::Int, has_header::Bool)

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
function generate_oscillation_probabilities(energies)

    cos_θ    = range(-1, 0, length = n_vectors)   # matches path count
    
    #energies, probs, paths=linkWithNeurthino()
    #export energies, probs, paths
    println("linkWithNeurthinoPREM works")
    load_first_interpolated_density_model()
    probs2, probs2anti = linkWithNeurthinoPREM_YD(cos_θ,energies)

    return probs2, probs2anti 

end

#PREMlineDensityElectron2D(fi, n_pts, iTime, detector,source, colorname, ax1, dR)
#import .Neurthino
#using Neurthino: Electron, Muon, Tau


#module Neutrino_Oscillations
#export generate_oscillation_probabilities

#using DIVAnd: DIVAnd_rectdom, DIVAndrun
#using Interpolations
#using DelimitedFiles
#using CairoMakie
#using DIVAnd
#using Colors
#using DelimitedFiles
#using JSON
#using FilePathsBase

###
###import .Neurthino
###using Neurthino: Electron, Muon, Tau

#ParamFile = "../config/testparam.csv" 
###include(joinpath(FLEXOPT_DIR, "src", "planet1D.jl"))
###using .planet1D
#using .DSM1D (?) Why is DSM1D relevant here?

###include(joinpath(FLEXOPT_DIR, "src","batchFiles", "batchUseful.jl"))
###include(joinpath(FLEXOPT_DIR, "src","batchFiles", "batchStagYY.jl"))
###include(joinpath(FLEXOPT_DIR, "src", "Neurthino", "NeurthinoRelated.jl"))
###include(joinpath(FLEXOPT_DIR, "src", "NeurthinoBack", "usefulFunctionsToPlot.jl"))
###include(joinpath(FLEXOPT_DIR, "src", "NeurthinoBack", "premFunctions.jl"))

# The following variables are taken from bare main, they should be passed as arguments in
# the functions within LinkWithNeurthinoPREM_YD


#### Cartesian grids and interpolation
###minX,maxX,nX = -6500e3, 6500e3, 521
###minY,maxY,nY = minX,maxX,nX
###minZ,maxZ,nZ = minX,maxX,nX
###dR = (maxX-minX)/(nX-1) # the interval in X, which we suppose to be the smallest grid interval

###correlationLength=(20e3,20e3,20e2) # not yet fully understood this for DIV interpolation

###epsilon2 =1.;

###nZ=1
###minZ=0.0
###maxZ=0.0
###tmpX=correlationLength[1]
###tmpY=correlationLength[2]
###correlationLength=(tmpX,tmpY)

###mask,(pm,pn),(xi,yi) = DIVAnd_rectdom(range(minX,stop=maxX,length=nX),
   ###                                     range(minY,stop=maxY,length=nY));
  
###iTime      = 2
###n_pts      = 100
###n_vectors  = 100      # ← changed

###zposition  = 2.5e3

###dir=joinpath(FLEXOPT_DIR, "op_old_full_mars_2025")
###rhoFiles=myListDir(dir; pattern=r"test_rho\d");
###compositionFiles=myListDir(dir; pattern=r"test_c\d");
###temperatureFiles=myListDir(dir; pattern=r"test_t\d");
###wtrFiles=myListDir(dir; pattern=r"test_wtr\d");

###rhoFiles = filter(f -> !occursin(r"/\._", f), rhoFiles) #if op_old_full_mars_2025



    #= Neurthino tests

    function creationPaths(n_vectors, zposition)

        densities_list, sections_list = vectorsFromDetector(n_vectors, zposition) 
        paths = Vector{Path}(undef, n_vectors)  

        for i in eachindex(paths)
            paths[i]= Path(densities_list[i],sections_list[i])
        end

        return paths
    end
    =#


#end 



  # path to flexOPT directory -> TO-DO: adapt this path to FlexOPT
    #include(joinpath(FLEXOPT_DIR, "src", "batchFiles",    "batchUseful.jl"))
    #include(joinpath(FLEXOPT_DIR, "src", "batchFiles",    "batchStagYY.jl"))
    #include(joinpath(FLEXOPT_DIR, "src", "Neurthino",     "NeurthinoRelated.jl"))
    #include(joinpath(FLEXOPT_DIR, "src", "NeurthinoBack", "usefulFunctionsToPlot.jl"))
    #include(joinpath(FLEXOPT_DIR, "src", "NeurthinoBack", "premFunctions.jl"))
    #dir=joinpath(FLEXOPT_DIR, "op_old_full_mars_2025")
    #rhoFiles=myListDir(dir; pattern=r"test_rho\d");
    #compositionFiles=myListDir(dir; pattern=r"test_c\d");
    #temperatureFiles=myListDir(dir; pattern=r"test_t\d");
    #wtrFiles=myListDir(dir; pattern=r"test_wtr\d");
    #rhoFiles = filter(f -> !occursin(r"/\._", f), rhoFiles)




       # number of paths
    n_paths   = 100   
    # number of points on each path
    n_pts     = 100
    # depth of the detector
    zposition = 2.5e3
dR = (maxX-minX)/(nX-1) # the smallest grid interval in X
#fieldname = "rho"
