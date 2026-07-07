const FLEXOPT_DIR = "/Users/igoos/Desktop/projects/flexOPT" # TO-DO: adapt this path to FlexOPT
include(joinpath(FLEXOPT_DIR, "src", "Neurthino.jl"))

Neurthino.setθ!(osc, 1=>2, 0.59)
Neurthino.setθ!(osc, 1=>3, 0.15)
Neurthino.setθ!(osc, 2=>3, 0.84)
Neurthino.setδ!(osc, 1=>3, 3.86)
Neurthino.setΔm²!(osc, 2=>3, -2.523e-3)
Neurthino.setΔm²!(osc, 1=>2, -7.39e-5)

"""
    set_oscillation_parameters(ordering=:normal)

Initializes a 3-flavor neutrino oscillation object 
using the latest global best-fit values from NuFIT v6.1.
Accepts either `:normal` or `:inverted` for the `ordering` keyword.
"""
function set_oscillation_parameters(; ordering::Symbol=:normal)
 
    # initialize a 3-flavor oscillation object
    osc = Neurthino.OscillationParameters(3)
    
    if ordering === :normal
        # --- NuFIT v6.1 Normal Ordering (Best Fit) ---
        # Mixing angles (converted to radians)
        Neurthino.setθ!(osc, 1=>2, 0.588)   # θ12 ≈ 33.68°
        Neurthino.setθ!(osc, 1=>3, 0.149)   # θ13 ≈ 8.56°
        Neurthino.setθ!(osc, 2=>3, 0.756)   # θ23 ≈ 43.3° (with SK-atm data)
        
        # CP-violating phase (radians)
        Neurthino.setδ!(osc, 1=>3, 3.700)   # δCP ≈ 212°
        
        # Mass-squared splittings (eV²)
        Neurthino.setΔm²!(osc, 1=>2,  7.49e-5)   # Δm²21
        Neurthino.setΔm²!(osc, 2=>3,  2.513e-3)  # Δm²31
        
    elseif ordering === :inverted
        # --- NuFIT v6.1 Inverted Ordering ---
        # Mixing angles (converted to radians)
        Neurthino.setθ!(osc, 1=>2, 0.588)   # θ12 ≈ 33.68°
        Neurthino.setθ!(osc, 1=>3, 0.150)   # θ13 ≈ 8.59°
        Neurthino.setθ!(osc, 2=>3, 0.836)   # θ23 ≈ 47.9° 
        
        # CP-violating phase (radians)
        Neurthino.setδ!(osc, 1=>3, 4.782)   # δCP ≈ 274°
        
        # Mass-squared splittings (eV²)
        Neurthino.setΔm²!(osc, 1=>2,  7.49e-5)   # Δm²21
        Neurthino.setΔm²!(osc, 2=>3, -2.484e-3)  # Δm²32 (Negative for IO)
        
    else
        error("Unknown mass ordering: :\$ordering. Use :normal or :inverted.")
    end
    
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
