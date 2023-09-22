function plot_gp_kernel_lattice_1d_lines(;z::Vector{Float64}=[0.,0.,0.,0.],o::Vector{Int64}=[2,2,2,2],γ::Vector{Float64}=[1.,1.,1.,1.],η::Vector{Float64}=[1.,1.,1.,1.],partial_x_order::Vector{Int64}=[0,1,0,1],partial_z_order::Vector{Int64}=[0,0,1,1],xmin::Float64=-.1,xmax::Float64=1.1,nxticks::Int64=257,linewidth::Float64=3.,backgroundcolor::Symbol=:white)
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

function plot_gp_kernel_lattice_1d_contour(;o::Vector{Int64}=[2,2,2,2],γ::Vector{Float64}=[1.,1.,1.,1.],η::Vector{Float64}=[1.,1.,1.,1.],partial_x_order::Vector{Int64}=[0,1,0,1],partial_z_order::Vector{Int64}=[0,0,1,1],xmin::Float64=-.1,xmax::Float64=1.1,nxticks::Int64=129,backgroundcolor::Symbol=:white)
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

function plot_gp_optimization(gp::GaussianProcessLatticeSeqB2;backgroundcolor::Symbol=:white)
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

function plot_gp_1s(gp::GaussianProcessLatticeSeqB2;f::Union{Nothing,Function}=nothing,partial_order::Vector{Int64}=[0],uncertainty::Float64=.05,xmin::Float64=-.1,xmax::Float64=1.1,nxticks::Int64=257,linewidth::Float64=3.,markersize::Float64=15.,backgroundcolor::Symbol=:white)
    n = length(partial_order)
    fig = CairoMakie.Figure(resolution=(800,n*500),backgroundcolor=backgroundcolor)
    xticks = Vector(xmin:(xmax-xmin)/(nxticks-1):xmax)
    beta = quantile(Normal(),1-uncertainty/2)
    if f!==nothing yticks = reshape(vcat([f([xticks[i]]) for i=1:nxticks]'...),nxticks,gp.r) end 
    for i=1:n
        ax = CairoMakie.Axis(fig[2*i,1],xlabel=L"$x$")
        CairoMakie.xlims!(ax,xmin,xmax)
        po = partial_order[i]
        idx = findfirst(x->x==po,gp.partial_orders[:,1])
        if (f!==nothing)&&(idx!==nothing) CairoMakie.lines!(ax,xticks,yticks[:,idx],color=JULIA4LOGOCOLORS[2],linewidth=linewidth,label=latexstring("\$f^{($po)}(x)\$")) end 
        if idx!==nothing CairoMakie.scatter!(ax,gp.x[:,1],gp.y[:,idx],markersize=markersize,color=JULIA4LOGOCOLORS[4],label=latexstring("\$(y^{($po)}_i)_{i=1}^{$(gp.n)}\$")) end 
        yhatticks = map(xtick->gp([xtick],[po]),xticks)
        stdhatticks = sqrt.(map(xtick->var_post(gp,[xtick],[po]),xticks))
        ci_low,ci_high = yhatticks.-beta*stdhatticks,yhatticks.+beta*stdhatticks
        CairoMakie.lines!(ax,xticks,yhatticks,color=JULIA4LOGOCOLORS[1],linewidth=linewidth,label=latexstring("\$m_n^{($po)}(x)\$"))
        CairoMakie.band!(ax,xticks,ci_low,ci_high,color=(JULIA4LOGOCOLORS[1],.25),label=latexstring("\$m_n^{($po)}(x) \\pm $(round(beta,digits=2)) \\; \\sigma_n^{($pi)}(x)\$"))
        CairoMakie.Legend(fig[2*i-1,1],ax,orientation=:horizontal,framevisible=false) 
    end 
    fig 
end 
