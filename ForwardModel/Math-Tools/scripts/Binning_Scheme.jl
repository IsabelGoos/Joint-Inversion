#module Binning_Scheme
#export set_binning_scheme

using CairoMakie


function set_energy_binning(Eν_min, Eν_max, nbins , scheme::Symbol =:linear ; return_edges::Bool  = false )

    if Eν_min <= 0
        error("Lower energy bound (Ev_min=$Eν_min) is not physical, must be > 0 ")
    end


    if scheme === :linear
        uEdges = range(Eν_min, Eν_max, nbins+1) # First element is the lower edge of the first bin 
     elseif scheme === :log10
        uEdges = logrange(Eν_min, Eν_max, nbins+1) # First element is the lower edge of the first bin 
     else
        error("$scheme is not known scheme. Options are linear or Log10, otherwise use custom_energy_binning") 
    end

    ΔEν =  diff(uEdges)
    Eν = 0.5 .* (uEdges[1:end-1] .+ uEdges[2:end]) # Bin centers

    return return_edges ? (Eν, ΔEν, uEdges) : (Eν, ΔEν)
end

function set_angular_binning(θ_min, θ_max, nbins , scheme::Symbol =:linear ; return_edges::Bool  = false )

    #TO-DO: theta must be in radias?
  
    if !( 0 <= θ_min < θ_max <= π )

        error("Bounds are not physical. For neutrino earth propagation analysis 0 <= θ_min < θ_max <= π (radians), " *
              "Perhaps you are using degrees instead of radians?")
    
    end

    #Linear bins in \theta variable
    if scheme === :linear #Checked :)

        uEdges = range(θ_min, θ_max, nbins+1) # First element is the lower edge of the first bin 
        
        cosEdges = reverse(cos.(uEdges)) # Recall that cos(\theta_max) < cos(\theta_min), then we need to re order bins
        
        Δcosθ = diff(cosEdges)

        θ_c = 0.5 .* (uEdges[1:end-1] .+ uEdges[2:end]) # zenith Bin centers?

        cosθ = reverse(cos.(θ_c))
        #cosθ = 0.5 .* (cosEdges[1:end-1] .+ cosEdges[2:end]) # cosine Bin centers?

        #Linear bins in cos\theta variable
    elseif scheme === :cos #Checked :)

            cos_min = cos(θ_max)
            cos_max = cos(θ_min)
            cosEdges = range(cos_min, cos_max, nbins+1) # First element is the lower edge of the first bin 
            Δcosθ = diff(cosEdges)
            cosθ = 0.5 .* (cosEdges[1:end-1] .+ cosEdges[2:end]) # Bin centers

            #Linear bins in log10(\theta) variable -> For core-crossing direction only!!!!

    elseif scheme === :log10angle
            tol = 0.01
            if (θ_min <= tol)
                @info "log10angle: θ_min ≤ tol, nudging lower edge to $tol rad to keep log spacing finite"
                θ_nudge =  tol 
                uEdges = logrange( θ_nudge, θ_max, nbins+1 )  
            else 
                uEdges = logrange( θ_min, θ_max, nbins+1)
            end

            cosEdges = reverse(cos.(uEdges))

            Δcosθ = diff(cosEdges)
            
            θ_c = 0.5 .* (uEdges[1:end-1] .+ uEdges[2:end]) # zenith Bin centers?

            cosθ = reverse(cos.(θ_c))

            #cosθ = 0.5 .* (cosEdges[1:end-1] .+ cosEdges[2:end]) # Bin centers

        #Linear bins in log10(cos\theta) variable -> For core-crossing direction only!!!!
    elseif scheme === :log10cos

                if !( π/2 <= θ_min < θ_max <= π )

                error("Bounds are not physical. The schemes is for earth propagation direction analysi so structly  π/2 <= θ_min < θ_max <= π (radians), " *
                    "Make sure tu use direction that correspond to neutrino propagation through the Earth   θ_min > π/2")
            
                end

            tol = 0.01 #Although it does not make sense to span the range [-1,0] tolerace for nudging, How close ftom 0 do we want   log10(0) \sim log10(tol)
            
            cos_min = cos(θ_max) # min vale is -1 (cos(π))

            cos_max = cos(θ_min) # max vale is 0 (cos(π/2))
            @info "log10cos: cos_min ≤ tol, nudging lower edge to $tol rad to keep log spacing finite"
            cosEdges = -reverse(logrange( max(abs(cos_max), tol),abs(cos_min), nbins+1 ))

            Δcosθ =  diff(cosEdges)

            cosθ = 0.5 .* (cosEdges[1:end-1] .+ cosEdges[2:end]) # Bin centers

    else
        error("$scheme is not known scheme. Options are linear or Log10, otherwise use custom_energy_binning") 
    end

    return return_edges ? (cosθ, Δcosθ, cosEdges) : (cosθ, Δcosθ)

end

# Binnig display

Eν_min = 1

Eν_max = 100

nbins = 10

Eν, ΔEν, uEdges = set_energy_binning(Eν_min,Eν_max,nbins,:linear; return_edges = true)


fig = Figure(size = (800, 400))
ax  = Axis(fig[1, 1], 
           xlabel = "Energy (GeV)",
           title  = "Energy bin edges, centers and widths",
           xscale = log10) 

# lower and upper edges as vertical lines
vlines!(ax, uEdges, color = :black, linestyle = :dash, label = "edges")

# centers as points
scatter!(ax, Eν, fill(1.0, length(Eν)),
         color = :red, markersize = 12, label = "centers")

# bin widths as horizontal spans
for i in 1:length(Eν)
    vspan!(ax, uEdges[i], uEdges[i+1],
           color = (:blue, 0.1 + 0.1 * (i % 2)))  # alternating shading
end

axislegend(ax, position = :lt)
ylims!(ax, 0, 2)        # y axis is meaningless here, just for visual space
hideydecorations!(ax)   # hide y axis since it carries no information
display(fig)

#%% Binning in zenith variable
θ_min, θ_max, nbins = 0, π, 10
cosθ, Δcosθ, cosEdges = set_angular_binning(θ_min, θ_max, nbins, :linear; return_edges = true)

angleEdges = reverse(acos.(cosEdges))   # ascending π/2 → π
θ_c        = reverse(acos.(cosθ))                # center angles (reverse not needed)
Δθ         = diff(angleEdges)           # correct angular widths (constant here)

fig = Figure(size = (800, 400))
ax  = Axis(fig[1, 1], xlabel = "θ (rad)",
           title = "Angular bin edges, centers and widths (scheme = :linear)")

vlines!(ax, angleEdges, color = :black, linestyle = :dash, label = "edges")
scatter!(ax, θ_c, fill(1.0, length(θ_c)), color = :red, markersize = 12, label = "centers")

for i in 1:length(θ_c)
    vspan!(ax, angleEdges[i], angleEdges[i+1], color = (:blue, 0.1 + 0.1 * (i % 2)))
end

axislegend(ax, position = :lt)
ylims!(ax, 0, 2)
hideydecorations!(ax)
display(fig)

#%% Cosine binning
θ_min = 0
θ_max = π
nbins = 10
cosθ, Δcosθ, cosEdges = set_angular_binning(θ_min, θ_max, nbins, :cos; return_edges = true)

fig = Figure(size = (800, 400))
ax  = Axis(fig[1, 1],
           xlabel = "cos θ",
           title  = "Angular bin edges, centers and widths (scheme = :cos)")
           # note: NO xscale = log10 — cosEdges are negative

# lower and upper edges as vertical lines
vlines!(ax, cosEdges, color = :black, linestyle = :dash, label = "edges")

# centers as points
scatter!(ax, cosθ, fill(1.0, length(cosθ)),
         color = :red, markersize = 12, label = "centers")

# bin widths as horizontal spans
for i in 1:length(cosθ)
    vspan!(ax, cosEdges[i], cosEdges[i+1],
           color = (:blue, 0.1 + 0.1 * (i % 2)))
end

axislegend(ax, position = :lt)
ylims!(ax, 0, 2)
hideydecorations!(ax)
display(fig)