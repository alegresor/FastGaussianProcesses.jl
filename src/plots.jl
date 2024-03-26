function plot_gp_fast_kernel_1s_lines(kernel_func::Function,x2::Vector{Float64},β1::Vector{Int64},β2::Vector{Int64},α::Vector{Int64},γ::Vector{Float64},η::Vector{Float64},xmin::Float64,xmax::Float64,nxticks::Int64,markersize::Float64,backgroundcolor::Symbol,figpath::Union{Nothing,String},px_per_unit::Int64)
    fig = CairoMakie.Figure(resolution=(800,400),backgroundcolor=backgroundcolor)
    ax = CairoMakie.Axis(fig[1,1], 
        xlabel = L"$x_1$",
        ylabel = L"$K^{(\beta_1,\beta_2)}(x_1,\; x_2 \; | \; \alpha \; \gamma, \; \eta)",
        backgroundcolor = backgroundcolor)
    xticks = Vector(xmin:(xmax-xmin)/nxticks:xmax)[1:end-1]
    n = length(γ); @assert n≤length(JULIA4LOGOCOLORS)
    for k=1:n
        kticks = [γ[k]*(β1[k]+β2[k] == 0 ? 1+η[k]*kernel_func(xticks[i],x2[k],β1[k],β2[k],α[k]) : η[k]*kernel_func(xticks[i],x2[k],β1[k],β2[k],α[k])) for i=1:nxticks]
        label = latexstring("\$x_2 = $(x2[k]), \\; \\beta_1 = $(β1[k]), \\; \\beta_2 = $(β2[k]), \\; \\alpha = $(α[k]), \\; \\gamma = $(γ[k]), \\; \\eta = $(η[k])\$")
        CairoMakie.scatter!(ax,xticks,kticks,color=JULIA4LOGOCOLORS[k],markersize=markersize,label=label)
    end 
    CairoMakie.Legend(fig[1,2],ax,framevisible=false)
    if figpath !== nothing CairoMakie.save(figpath,fig,px_per_unit=px_per_unit) end 
    fig
end
plot_gp_kernel_latticeseqb2_1s_lines(;x2::Vector{Float64}=[0.,0.,0.,0.],β1::Vector{Int64}=[0,1,0,1],β2::Vector{Int64}=[0,0,1,1],α::Vector{Int64}=[4,4,4,4],γ::Vector{Float64}=[1.,1.,1.,1.],η::Vector{Float64}=[1.,1.,1.,1.],xmin::Float64=-.1,xmax::Float64=1.1,nxticks::Int64=1024,markersize::Float64=8.,backgroundcolor::Symbol=:white,figpath::Union{Nothing,String}=nothing,px_per_unit::Int64=4) = plot_gp_fast_kernel_1s_lines(kernel_shiftinvar_s1,x2,β1,β2,α,γ,η,xmin,xmax,nxticks,markersize,backgroundcolor,figpath,px_per_unit)
plot_gp_kernel_digitalseqb2g_1s_lines(;x2::Vector{Float64}=[0.,0.,0.,0.],β1::Vector{Int64}=[0,1,0,1],β2::Vector{Int64}=[0,0,1,1],α::Vector{Int64}=[4,4,4,4],γ::Vector{Float64}=[1.,1.,1.,1.],η::Vector{Float64}=[1.,1.,1.,1.],xmin::Float64=0.,xmax::Float64=1.,nxticks::Int64=1024,markersize::Float64=8.,backgroundcolor::Symbol=:white,figpath::Union{Nothing,String}=nothing,px_per_unit::Int64=4) = plot_gp_fast_kernel_1s_lines(kernel_digshiftinvar_s1,x2,β1,β2,α,γ,η,xmin,xmax,nxticks,markersize,backgroundcolor,figpath,px_per_unit)
function plot_gp_kernel_rbf_1s_lines(;x2::Vector{Float64}=[0.,0.,0.,0.],β1::Vector{Int64}=[0,1,0,1],β2::Vector{Int64}=[0,0,1,1],γ::Vector{Float64}=[1.,1.,1.,1.],η::Vector{Float64}=[1.,1.,1.,1.],xmin::Float64=-3.5,xmax::Float64=3.5,nxticks::Int64=1024,markersize::Float64=8.,backgroundcolor::Symbol=:white,figpath::Union{Nothing,String}=nothing,px_per_unit::Int64=4)
    fig = CairoMakie.Figure(resolution=(800,400),backgroundcolor=backgroundcolor)
    ax = CairoMakie.Axis(fig[1,1], 
        xlabel = L"$x_1$",
        ylabel = L"$K^{(\beta_1,\beta_2)}(x_1,\; x_2 \; | \; \gamma, \; \eta)",
        backgroundcolor = backgroundcolor)
    xticks = Vector(xmin:(xmax-xmin)/nxticks:xmax)[1:end-1]
    n = length(γ); @assert n≤length(JULIA4LOGOCOLORS)
    for k=1:n
        kticks = [rbf_kernel([xticks[i]],[x2[k]],[β1[k]],[β2[k]],γ[k],[η[k]]) for i=1:nxticks]
        label = latexstring("\$x_2 = $(x2[k]), \\; \\beta_1 = $(β1[k]), \\; \\beta_2 = $(β2[k]), \\; \\gamma = $(γ[k]), \\; \\eta = $(η[k])\$")
        CairoMakie.scatter!(ax,xticks,kticks,color=JULIA4LOGOCOLORS[k],markersize=markersize,label=label)
    end 
    CairoMakie.Legend(fig[1,2],ax,framevisible=false)
    if figpath !== nothing CairoMakie.save(figpath,fig,px_per_unit=px_per_unit) end 
    fig
end 

function plot_gp_kernel_1s_contsurfs(kernel_func::Function,β1::Vector{Int64},β2::Vector{Int64},α::Vector{Int64},γ::Vector{Float64},η::Vector{Float64},xmin::Float64,xmax::Float64,nxticks::Int64,backgroundcolor::Symbol,figpath::Union{Nothing,String},px_per_unit::Int64)
    n = length(γ)
    fig = CairoMakie.Figure(resolution=(800,310*n),backgroundcolor=backgroundcolor)
    x1ticks = Vector(xmin:(xmax-xmin)/nxticks:xmax)[1:end-1]
    x2ticks = copy(x1ticks)
    for k=1:n 
        karr = [γ[k]*(β1[k]+β2[k] == 0 ? 1+η[k]*kernel_func(x1ticks[i],x2ticks[j],β1[k],β2[k],α[k]) : η[k]*kernel_func(x1ticks[i],x2ticks[j],β1[k],β2[k],α[k])) for i=1:nxticks,j=1:nxticks]
        kmin,kmax = minimum(karr),maximum(karr)
        ax = CairoMakie.Axis3(fig[k,1],
            xlabel = L"$x_1$",
            ylabel = L"$x_2$",
            zlabel = L"$K^{(\beta_1,\beta_2)}(x_1,\; x_2 \; |  \; \alpha, \; \gamma, \; \eta)",
            title = latexstring("\$\\beta_1 = $(β1[k]), \\; \\beta_2 = $(β2[k]), \\; \\alpha = $(α[k]), \\; \\gamma = $(γ[k]), \\; \\eta = $(η[k])\$"),
            backgroundcolor = backgroundcolor)
        CairoMakie.surface!(ax,x1ticks,x2ticks,karr,colormap=:julia_colorscheme,colorrange=(kmin,kmax))
        CairoMakie.xlims!(ax,xmin,xmax); CairoMakie.ylims!(ax,xmin,xmax); CairoMakie.zlims!(kmin,kmax)
        ax = CairoMakie.Axis(fig[k,2],
            xlabel = L"$x_1$",
            ylabel = L"$x_2$",
            aspect = 1,
            backgroundcolor = backgroundcolor)
        CairoMakie.heatmap!(ax,x1ticks,x2ticks,karr,colormap=:julia_colorscheme,colorrange=(kmin,kmax))
        CairoMakie.xlims!(ax,xmin,xmax); CairoMakie.ylims!(ax,xmin,xmax);
    end 
    if figpath !== nothing CairoMakie.save(figpath,fig,px_per_unit=px_per_unit) end 
    fig
end
plot_gp_kernel_latticeseqb2_1s_contsurfs(;β1::Vector{Int64}=[0,1,0,1],β2::Vector{Int64}=[0,0,1,1],α::Vector{Int64}=[2,2,2,2],γ::Vector{Float64}=[1.,1.,1.,1.],η::Vector{Float64}=[1.,1.,1.,1.],xmin::Float64=-.1,xmax::Float64=1.1,nxticks::Int64=128,backgroundcolor::Symbol=:white,figpath::Union{Nothing,String}=nothing,px_per_unit::Int64=4) = plot_gp_kernel_1s_contsurfs(kernel_shiftinvar_s1,β1,β2,α,γ,η,xmin,xmax,nxticks,backgroundcolor,figpath,px_per_unit)
plot_gp_kernel_digitalseqb2g_1s_contsurfs(;β1::Vector{Int64}=[0,1,1],β2::Vector{Int64}=[0,0,1],α::Vector{Int64}=[4,4,4],γ::Vector{Float64}=[1.,1.,1.],η::Vector{Float64}=[1.,1.,1.],xmin::Float64=0.,xmax::Float64=1.,nxticks::Int64=128,backgroundcolor::Symbol=:white,figpath::Union{Nothing,String}=nothing,px_per_unit::Int64=4) = plot_gp_kernel_1s_contsurfs(kernel_digshiftinvar_s1,β1,β2,α,γ,η,xmin,xmax,nxticks,backgroundcolor,figpath,px_per_unit)
function plot_gp_kernel_rbf_1s_contsurfs(;β1::Vector{Int64}=[0,1,0,1],β2::Vector{Int64}=[0,0,1,1],γ::Vector{Float64}=[1.,1.,1.,1.],η::Vector{Float64}=[1.,1.,1.,1.],xmin::Float64=-3.5,xmax::Float64=3.5,nxticks::Int64=128,backgroundcolor::Symbol=:white,figpath::Union{Nothing,String}=nothing,px_per_unit::Int64=4)
    n = length(γ)
    fig = CairoMakie.Figure(resolution=(800,310*n),backgroundcolor=backgroundcolor)
    x1ticks = Vector(xmin:(xmax-xmin)/nxticks:xmax)[1:end-1]
    x2ticks = copy(x1ticks)
    for k=1:n 
        karr = [rbf_kernel([x1ticks[i]],[x2ticks[j]],[β1[k]],[β2[k]],γ[k],[η[k]]) for i=1:nxticks,j=1:nxticks]
        kmin,kmax = minimum(karr),maximum(karr)
        ax = CairoMakie.Axis3(fig[k,1],
            xlabel = L"$x_1$",
            ylabel = L"$x_2$",
            zlabel = L"$K^{(\beta_1,\beta_2)}(x_1,\; x_2 \; | \; \gamma, \; \eta)",
            title = latexstring("\$\\beta_1 = $(β1[k]), \\; \\beta_2 = $(β2[k]), \\; \\gamma = $(γ[k]), \\; \\eta = $(η[k])\$"),
            backgroundcolor = backgroundcolor)
        CairoMakie.surface!(ax,x1ticks,x2ticks,karr,colormap=:julia_colorscheme,colorrange=(kmin,kmax))
        CairoMakie.xlims!(ax,xmin,xmax); CairoMakie.ylims!(ax,xmin,xmax); CairoMakie.zlims!(kmin,kmax)
        ax = CairoMakie.Axis(fig[k,2],
            xlabel = L"$x_1$",
            ylabel = L"$x_2$",
            aspect = 1,
            backgroundcolor = backgroundcolor)
        CairoMakie.heatmap!(ax,x1ticks,x2ticks,karr,colormap=:julia_colorscheme,colorrange=(kmin,kmax))
        CairoMakie.xlims!(ax,xmin,xmax); CairoMakie.ylims!(ax,xmin,xmax);
    end 
    if figpath !== nothing CairoMakie.save(figpath,fig,px_per_unit=px_per_unit) end 
    fig
end 

function plot_gp_optimization(gp::Union{FastGaussianProcess,GaussianProcessRBF};backgroundcolor::Symbol=:white,figpath::Union{Nothing,String}=nothing,px_per_unit::Int64=4)
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
    for l=1:gp.n_β CairoMakie.lines!(axζ,xrange,gp.ζs[:,l],color=(JULIA4LOGOCOLORS[4],(gp.n_β-l+1)/gp.n_β),linewidth=3) end 
    CairoMakie.xlims!(axloss,0,noptsp1-1); CairoMakie.xlims!(axζ,0,noptsp1-1); CairoMakie.xlims!(axη,0,noptsp1-1); CairoMakie.xlims!(axζ,0,noptsp1-1)
    if figpath !== nothing CairoMakie.save(figpath,fig,px_per_unit=px_per_unit) end 
    return fig
end

function plot_gp_1s(gp::Union{FastGaussianProcess,GaussianProcessRBF};f::Union{Nothing,Function}=nothing,β::Vector{Int64}=[0],uncertainty::Float64=.05,xmin::Float64=0.,xmax::Float64=1.,nxticks::Int64=1024,markersize::Float64=16.,backgroundcolor::Symbol=:white,figpath::Union{Nothing,String}=nothing,px_per_unit::Int64=4)
    @assert gp.s==1 
    n = length(β)
    fig = CairoMakie.Figure(resolution=(800,n*500),backgroundcolor=backgroundcolor)
    xticks = Vector(xmin:(xmax-xmin)/nxticks:xmax)[1:end-1]
    q = quantile(Normal(),1-uncertainty/2)
    if f!==nothing yticks = reshape(vcat([f([xticks[i]]) for i=1:nxticks]'...),nxticks,gp.n_β) end 
    for i=1:n
        ax = CairoMakie.Axis(fig[2*i,1],xlabel=L"$x$")
        CairoMakie.xlims!(ax,xmin,xmax)
        po = β[i]
        idx = findfirst(x->x==po,gp.β[:,1])
        if (f!==nothing)&&(idx!==nothing) CairoMakie.scatter!(ax,xticks,yticks[:,idx],color=JULIA4LOGOCOLORS[2],markersize=markersize/2,label=latexstring("\$f^{($po)}(x)\$")) end 
        yhatticks = map(xtick->gp([xtick];β=[po]),xticks)
        CairoMakie.scatter!(ax,xticks,yhatticks,color=JULIA4LOGOCOLORS[1],markersize=markersize/2,label=latexstring("\$m_n^{($po)}(x)\$"))
        stdhatticks = sqrt.(map(xtick->var_post(gp,[xtick];β=[po]),xticks))
        ci_low,ci_high = yhatticks.-q*stdhatticks,yhatticks.+q*stdhatticks
        CairoMakie.band!(ax,xticks,ci_low,ci_high,color=(JULIA4LOGOCOLORS[1],.25),label=latexstring("\$m_n^{($po)}(x) \\pm $(round(q,digits=2)) \\; \\sigma_n^{($po)}(x)\$"))
        if idx!==nothing CairoMakie.scatter!(ax,gp.x[:,1],gp.y[:,idx],markersize=markersize,color=:black,label=latexstring("\$(y^{($po)}_i)_{i=1}^{$(gp.n)}\$")) end 
        CairoMakie.Legend(fig[2*i-1,1],ax,orientation=:horizontal,framevisible=false) 
    end
    if figpath !== nothing CairoMakie.save(figpath,fig,px_per_unit=px_per_unit) end 
    fig 
end

function plot_gp_2s(gp::Union{FastGaussianProcess,GaussianProcessRBF};f::Union{Nothing,Function}=nothing,β::Matrix{Int64}=[0 0;],xmin::Float64=0.,xmax::Float64=1.,nxticks::Int64=32,markersize::Float64=4.,backgroundcolor::Symbol=:white,figpath::Union{Nothing,String}=nothing,px_per_unit::Int64=4)
    @assert gp.s==2
    n = size(β,1)
    cols = f===nothing ? 1 : 2
    fig = CairoMakie.Figure(resolution=(cols*2*400,n*300),backgroundcolor=backgroundcolor)
    xticks = Vector(xmin:(xmax-xmin)/nxticks:xmax)[1:end-1]
    if f!==nothing yticks = reshape([vcat(f([xticks[i],xticks[j]])) for i=1:nxticks,j=1:nxticks],nxticks,nxticks) end 
    for i=1:n
        po = β[i,:]; po1,po2 = po[1],po[2]
        idx = findfirst(x->x==1,all(gp.β.==reshape(po,1,2),dims=2))
        if (f!==nothing)&&(idx!==nothing)
            ymesh = [yticks[i,j][idx] for i=1:nxticks,j=1:nxticks]
            ax = CairoMakie.Axis(fig[i,1],xlabel=L"$x_1$",ylabel=L"$x_2$",aspect=1,title=latexstring("\$f^{($po1,$po2)}(x)\$")); CairoMakie.xlims!(ax,xmin,xmax); CairoMakie.ylims!(ax,xmin,xmax)
            CairoMakie.heatmap!(ax,xticks,xticks,ymesh,colormap=:julia_colorscheme) 
            CairoMakie.scatter!(ax,gp.x[:,1],gp.x[:,2],markersize=markersize,color=:black,label=latexstring("\$(y^{($po1,$po2)}_i)_{i=1}^{$(gp.n)}\$")) 
            ax = CairoMakie.Axis3(fig[i,2],xlabel=L"$x_1$",ylabel=L"$x_2$",zlabel="",title=latexstring("\$f^{($po1,$po2)}(x)\$")); CairoMakie.xlims!(ax,xmin,xmax); CairoMakie.ylims!(ax,xmin,xmax)
            CairoMakie.surface!(ax,xticks,xticks,ymesh,colormap=:julia_colorscheme) 
            CairoMakie.scatter!(ax,gp.x[:,1],gp.x[:,2],gp.y[:,idx],markersize=markersize,color=:black,label=latexstring("\$(y^{($po1,$po2)}_i)_{i=1}^{$(gp.n)}\$")) 
        end 
        yhatticks = [gp([xticks[i],xticks[j]];β=po) for i=1:nxticks,j=1:nxticks]
        ax = CairoMakie.Axis3(fig[i,f===nothing ? 1 : 3],xlabel=L"$x_1$",ylabel=L"$x_2$",zlabel="",title=latexstring("\$m_n^{($po1,$po2)}(x)\$")); CairoMakie.xlims!(ax,xmin,xmax); CairoMakie.ylims!(ax,xmin,xmax)
        CairoMakie.surface!(ax,xticks,xticks,yhatticks,colormap=:julia_colorscheme)
        if idx!==nothing CairoMakie.scatter!(ax,gp.x[:,1],gp.x[:,2],gp.y[:,idx],markersize=markersize,color=:black,label=latexstring("\$(y^{($po1,$po2)}_i)_{i=1}^{$(gp.n)}\$")) end 
        ax = CairoMakie.Axis(fig[i,f===nothing ? 2 : 4],xlabel=L"$x_1$",ylabel=L"$x_2$",aspect=1,title=latexstring("\$m_n^{($po1,$po2)}(x)\$")); CairoMakie.xlims!(ax,xmin,xmax); CairoMakie.ylims!(ax,xmin,xmax)
        CairoMakie.heatmap!(ax,xticks,xticks,yhatticks,colormap=:julia_colorscheme)
        if idx!==nothing CairoMakie.scatter!(ax,gp.x[:,1],gp.x[:,2],markersize=markersize,color=:black,label=latexstring("\$(y^{($po1,$po2)}_i)_{i=1}^{$(gp.n)}\$")) end 
    end
    if figpath !== nothing CairoMakie.save(figpath,fig,px_per_unit=px_per_unit) end 
    fig
end 
