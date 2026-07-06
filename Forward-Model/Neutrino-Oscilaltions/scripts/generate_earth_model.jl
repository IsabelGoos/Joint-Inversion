module Earth_Model
export read_earth_model

using DIVAnd

function read_earth_model()

    # cartesian grid in 2D
    minX, maxX, nX = -6500e3, 6500e3, 521 # m, m, number
    minY, maxY, nY = minX, maxX, nX
    dR = (maxX-minX)/(nX-1) # the smallest grid interval in X

    # interpolation
    # This tells the function how far apart two data points can be before they stop influencing each other.
    correlationLength = (20e3, 20e3) # m
    # This tells the function how much to trust your sensors. High trust -> low value
    epsilon2 = 1.0 
    mask, (pm, pn), (xi, yi) = DIVAnd_rectdom(range(minX, maxX, length=nX), range(minY,maxY, length=nY));
  
    # number of files
    iTime     = 2
    # number of paths
    n_paths   = 100   
    # number of points on each path
    n_pts     = 100
    # depth of the detector
    zposition = 2.5e3

    # path to flexOPT directory -> TO-DO: adapt this path to FlexOPT
    FLEXOPT_DIR = "/Users/igoos/Desktop/projects/flexOPT" 

end

end