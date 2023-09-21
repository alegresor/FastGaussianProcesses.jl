function plot_gp_optimization(gp::GaussianProcessLatticeSeqB2,backgroundcolor::Symbol)
    noptsp1 = length(gp.losses)
    @assert noptsp1>1
    xrange = [k for k=0:noptsp1-1]
    fig = CairoMakie.Figure(resolution=(800,600),sharex=true,backgroundcolor=backgroundcolor)
    axloss = CairoMakie.Axis(fig[1,1],
        ylabel = L"$\text{Loss } L(\theta | Y)$",
        backgroundcolor = backgroundcolor)
    CairoMakie.lines!(axloss,xrange,gp.losses,color=JULIA4LOGOCOLORS[1],linewidth=3)
    CairoMakie.xlims!(axloss,0,noptsp1-1)
    axγ = CairoMakie.Axis(fig[1,2],
        ylabel = L"$\gamma$",
        yscale = log10,
        backgroundcolor = backgroundcolor)
    CairoMakie.lines!(axγ,xrange,gp.γs,color=JULIA4LOGOCOLORS[2],linewidth=3)
    CairoMakie.xlims!(axγ,0,noptsp1-1)
    axη = CairoMakie.Axis(fig[2,1],
        xlabel = "step",
        ylabel = L"$\eta$",
        yscale = log10,
        backgroundcolor=backgroundcolor)
    for j=1:gp.s CairoMakie.lines!(axη,xrange,gp.ηs[:,j],color=(JULIA4LOGOCOLORS[3],(gp.s-j+1)/gp.s),linewidth=3) end 
    CairoMakie.xlims!(axη,0,noptsp1-1)
    axζ = CairoMakie.Axis(fig[2,2],
        xlabel = "step",
        ylabel = L"$\zeta$",
        yscale = log10,
        backgroundcolor = backgroundcolor)
    for l=1:gp.r CairoMakie.lines!(axζ,xrange,gp.ζs[:,l],color=(JULIA4LOGOCOLORS[4],(gp.r-l+1)/gp.r),linewidth=3) end 
    CairoMakie.xlims!(axζ,0,noptsp1-1)
    CairoMakie.linkxaxes!(axloss,axη)
    CairoMakie.linkxaxes!(axγ,axζ)
    return fig
end

function plot_gp_kernel_lattice_1d_lines(z::Vector{Float64},o::Vector{Int64},γ::Vector{Float64},η::Vector{Float64},partial_x_order::Vector{Int64},partial_z_order::Vector{Int64},xmin::Float64,xmax::Float64,nxticks::Int64,linewidth::Float64,backgroundcolor::Symbol)
    fig = CairoMakie.Figure(resolution=(800,400),backgroundcolor=backgroundcolor)
    ax = CairoMakie.Axis(fig[1,1],
        xlabel = L"$x$",
        ylabel = L"$k^{(\alpha,\beta)}(x,\; z \; |  \; o, \; \gamma, \; \eta)",
        backgroundcolor = backgroundcolor)
    xticks = Vector(xmin:(xmax-xmin)/(nxticks-1):xmax)
    n = length(o); @assert n≤length(JULIA4LOGOCOLORS)
    for k=1:n
        kticks = map(xtick->kernel_lattice([xtick],[z[k]],o[k],γ[k],[η[k]],[partial_x_order[k]],[partial_z_order[k]],1),xticks)
        label = latexstring("\$\\alpha = $(partial_x_order[k]), \\; \\beta = $(partial_z_order[k]), \\; z = $(z[k]), \\; o = $(o[k]), \\; \\gamma = $(γ[k]), \\; \\eta = $(η[k])\$")
        CairoMakie.lines!(ax,xticks,kticks,color=JULIA4LOGOCOLORS[k],linewidth=linewidth,label=label)
    end 
    CairoMakie.Legend(fig[1,2],ax,framevisible=false)
    fig
end

function plot_gp_kernel_lattice_1d_contour(o::Vector{Int64},γ::Vector{Float64},η::Vector{Float64},partial_x_order::Vector{Int64},partial_z_order::Vector{Int64},xmin::Float64,xmax::Float64,nxticks::Int64,backgroundcolor::Symbol)
    n = length(o)
    fig = CairoMakie.Figure(resolution=(800,310*n),backgroundcolor=backgroundcolor)
    x1ticks = Vector(xmin:(xmax-xmin)/(nxticks-1):xmax)
    x2ticks = copy(x1ticks)
    karr = [kernel_lattice([x1ticks[i]],[x2ticks[j]],o[k],γ[k],[η[k]],[partial_x_order[k]],[partial_z_order[k]],1) for k=1:n,i=1:nxticks,j=1:nxticks]
    for k=1:n 
        kmin,kmax = minimum(karr[k,:,:]),maximum(karr[k,:,:])
        ax = CairoMakie.Axis3(fig[k,1],
            xlabel = L"$x_1$",
            ylabel = L"$x_2$",
            zlabel = L"$k^{(\alpha,\beta)}(x_1,\; x_2 \; |  \; o, \; \gamma, \; \eta)",
            title = latexstring("\$\\alpha = $(partial_x_order[k]), \\; \\beta = $(partial_z_order[k]), \\; o = $(o[k]), \\; \\gamma = $(γ[k]), \\; \\eta = $(η[k])\$"),
            backgroundcolor = backgroundcolor)
        CairoMakie.surface!(ax,x1ticks,x2ticks,karr[k,:,:],colormap=:julia_colorscheme,colorrange=(kmin,kmax))
        CairoMakie.xlims!(ax,xmin,xmax); CairoMakie.ylims!(ax,xmin,xmax); CairoMakie.zlims!(kmin,kmax)
        ax = CairoMakie.Axis(fig[k,2],
            xlabel = L"$x_1$",
            ylabel = L"$x_2$",
            aspect = 1,
            backgroundcolor = backgroundcolor)
        CairoMakie.heatmap!(ax,x1ticks,x2ticks,karr[k,:,:],colormap=:julia_colorscheme,colorrange=(kmin,kmax))
        CairoMakie.xlims!(ax,xmin,xmax); CairoMakie.ylims!(ax,xmin,xmax);
    end 
    fig
end 