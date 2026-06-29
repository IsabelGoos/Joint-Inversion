using JSON

module Neutrino_Cross_Sections

function read_neutrino_cross_sections(filename::String)

    # construct path
    data_dir  = joinpath(@__DIR__, "..", "data")
    data_path = joinpath(data_dir, "$(filename)")

    # read the data
    if endswith(filename, ".json")
        data = JSON.parsefile(data_path)
    end

    struct CubicSpline_Info
        # knots
        k::Vector{Float64}
        # coefficients
        c::Matrix{Float64}
    end

    cross_section_splines = Dict{String, CubicSpline_Info}()
    for (key, val) in data
        k = Float64.(val["k"])
        c_list = val["c"]                
        c = Float64.(hcat(c_list...))'   
        cross_section_splines[key] = CubicSpline_Info(k, c)
    end

    nue_CC  = cross_section_splines["nue_CC_logE"]
    num_CC  = cross_section_splines["num_CC_logE"]
    nut_CC  = cross_section_splines["nut_CC_logE"]
    anue_CC = cross_section_splines["anue_CC_logE"]
    anum_CC = cross_section_splines["anum_CC_logE"]
    anut_CC = cross_section_splines["anut_CC_logE"]

    return nue_CC, num_CC, nut_CC, anue_CC, anum_CC, anut_CC

end

end


"""

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

"""