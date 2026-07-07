
import .Neurthino
using Neurthino: Electron, Muon, Tau


#module Neutrino_Oscillations
#export generate_oscillation_probabilities

using DIVAnd: DIVAnd_rectdom, DIVAndrun
#using Interpolations
#using DelimitedFiles
#using CairoMakie
#using DIVAnd
#using Colors
#using DelimitedFiles
#using JSON
#using FilePathsBase

# TO-DO: remove this path to FlexOPT 
###const FLEXOPT_DIR = "/Users/igoos/Desktop/projects/flexOPT" #Path to flexOPT directory
###include(joinpath(FLEXOPT_DIR, "src", "Neurthino.jl"))
###import .Neurthino
###using Neurthino: Electron, Muon, Tau

ParamFile = "../config/testparam.csv" 
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

    probs2, probs2anti = linkWithNeurthinoPREM_YD(cos_θ,energies)

    return probs2, probs2anti 

end

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