module Earth_Model
export read_earth_model

const FLEXOPT_DIR = "/Users/igoos/Desktop/projects/flexOPT" # TO-DO: adapt this path to FlexOPT
include(joinpath(FLEXOPT_DIR, "src", "batchFiles", "batchStagYY.jl"))
include(joinpath(FLEXOPT_DIR, "src", "NeurthinoBack", "premFunctions.jl"))
ParamFile = "../config/testparam.csv" 
include(joinpath(FLEXOPT_DIR, "src", "planet1D.jl"))
using DIVAnd

"""
    load_first_interpolated_density_model(; iTime=3)

Loads a 2D staggered-grid mantle density dataset, extends with PREM core density values,
and interpolates it onto a uniform 13,000 km x 13,000 km Cartesian grid.
Returns a 2D matrix of the interpolated density field.
"""
function load_first_interpolated_density_model(; iTime::Int=3)

    # define the 2D cartesian grid
    # create a 13,000 km square grid centered at (0,0) with 521 resolution nodes
    minX, maxX, nX = -6500e3, 6500e3, 521 # m, m, number
    minY, maxY, nY = minX, maxX, nX

    # calculate grid spacing resolution (distance between nodes)
    dR = (maxX-minX)/(nX-1) 

    # initialize the spatial domain metrics
    mask, (pm, pn), (xi, yi) = DIVAnd_rectdom(range(minX, maxX, length=nX), range(minY, maxY, length=nY))

    # read raw Staggered-Grid density data
    data_dir = joinpath(FLEXOPT_DIR, "op_old_full_mars_2025")
    rho_Files = myListDir(data_dir; pattern=r"test_rho\d")
    raw_field, Xnode, Ynode, _ = readStagYYFiles(rho_Files[iTime])

    # extend adding the core density from PREM since the rhoFiles don't have this
    extendWithρ!(raw_field, Xnode, Ynode, dR, iCheckCoreModel=false)

    # variational interpolation
    # This tells the function how far apart two data points can be before they stop influencing each other.
    correlationLength = (20e3, 20e3) # m
    # This tells the function how much to trust your sensors. High trust -> low value
    epsilon2 = 1.0 
    interpolated_density, _ = DIVAndrun(mask, (pm, pn), (xi, yi), (Xnode, Ynode), raw_field, correlationLength, epsilon2)

    return interpolated_density

end

end