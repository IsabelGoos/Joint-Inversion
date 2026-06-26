


# Esto es para leer el file, hay que acomodar el path
# Get cross sections
#data = JSON.parsefile("cross_sections_water.json")




"""
struct CubicSplineSciPy
    k::Vector{Float64}
    c::Matrix{Float64}
end
function evaluate(s::CubicSplineSciPy, x::Float64)
    k = s.k
    c = s.c
    # find interval
    i = searchsortedlast(k, x)
    if i >= length(k)
        i = length(k) - 1
    elseif i < 1
        i = 1
    end
    dx = x - k[i]
    return ((c[1,i]*dx + c[2,i])*dx + c[3,i])*dx + c[4,i]
end

cross_section_splines = Dict{String, CubicSplineSciPy}()
for (key, val) in data
    k = Float64.(val["k"])
    c_list = val["c"]                
    c = Float64.(hcat(c_list...))'   
    cross_section_splines[key] = CubicSplineSciPy(k, c)
end
"""

"""
nue_CC_logE = cross_section_splines["nue_CC_logE"]
num_CC_logE = cross_section_splines["num_CC_logE"]
nut_CC_logE = cross_section_splines["nut_CC_logE"]
anue_CC_logE = cross_section_splines["anue_CC_logE"]
anum_CC_logE = cross_section_splines["anum_CC_logE"]
anut_CC_logE = cross_section_splines["anut_CC_logE"]
"""