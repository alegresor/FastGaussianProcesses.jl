function plot_gp_kernel_1s_lines(kernel_func::Function,x2::Vector{Float64},β1::Vector{Int64},β2::Vector{Int64},α::Vector{Int64},γ::Vector{Float64},η::Vector{Float64},xmin::Float64,xmax::Float64,nxticks::Int64,markersize::Float64,backgroundcolor::Symbol)
    fig = CairoMakie.Figure(resolution=(800,400),backgroundcolor=backgroundcolor)
    ax = CairoMakie.Axis(fig[1,1], 
        xlabel = L"$x_1$",
        ylabel = L"$K^{(\beta_1,\beta_2)}(x_1,\; x_2 \; | \; \alpha \; \gamma, \; \eta)",
        backgroundcolor = backgroundcolor)
    xticks = Vector(xmin:(xmax-xmin)/nxticks:xmax)[1:end-1]
    n = length(α); @assert n≤length(JULIA4LOGOCOLORS)
    for k=1:n
        kticks = map(xtick->kernel_func([xtick],[x2[k]],[β1[k]],[β2[k]],α[k],γ[k],[η[k]],1),xticks)
        label = latexstring("\$x_2 = $(x2[k]), \\; \\beta_1 = $(β1[k]), \\; \\beta_2 = $(β2[k]), \\; \\alpha = $(α[k]), \\; \\gamma = $(γ[k]), \\; \\eta = $(η[k])\$")
        CairoMakie.scatter!(ax,xticks,kticks,color=JULIA4LOGOCOLORS[k],markersize=markersize,label=label)
    end 
    CairoMakie.Legend(fig[1,2],ax,framevisible=false)
    fig
end
plot_gp_kernel_latticeseqb2_1s_lines(;x2::Vector{Float64}=[0.,0.,0.,0.],β1::Vector{Int64}=[0,1,0,1],β2::Vector{Int64}=[0,0,1,1],α::Vector{Int64}=[2,2,2,2],γ::Vector{Float64}=[1.,1.,1.,1.],η::Vector{Float64}=[1.,1.,1.,1.],xmin::Float64=-.1,xmax::Float64=1.1,nxticks::Int64=1024,markersize::Float64=8.,backgroundcolor::Symbol=:white) = plot_gp_kernel_1s_lines(kernel_shiftinvar,x2,β1,β2,α,γ,η,xmin,xmax,nxticks,markersize,backgroundcolor)
plot_gp_kernel_digitalseqb2g_1s_lines(;x2::Vector{Float64}=[0.,0.,0.],β1::Vector{Int64}=[0,1,1],β2::Vector{Int64}=[0,0,1],α::Vector{Int64}=[4,4,4],γ::Vector{Float64}=[1.,1.,1.],η::Vector{Float64}=[1.,1.,1.],xmin::Float64=0.,xmax::Float64=1.,nxticks::Int64=1024,markersize::Float64=8.,backgroundcolor::Symbol=:white) = plot_gp_kernel_1s_lines(kernel_digshiftinvar,x2,β1,β2,α,γ,η,xmin,xmax,nxticks,markersize,backgroundcolor)

function plot_gp_kernel_1s_contsurfs(kernel_func::Function,β1::Vector{Int64},β2::Vector{Int64},α::Vector{Int64},γ::Vector{Float64},η::Vector{Float64},xmin::Float64,xmax::Float64,nxticks::Int64,backgroundcolor::Symbol)
    n = length(α)
    fig = CairoMakie.Figure(resolution=(800,310*n),backgroundcolor=backgroundcolor)
    x1ticks = Vector(xmin:(xmax-xmin)/nxticks:xmax)[1:end-1]
    x2ticks = copy(x1ticks)
    karr = [kernel_func([x1ticks[i]],[x2ticks[j]],[β1[k]],[β2[k]],α[k],γ[k],[η[k]],1) for k=1:n,i=1:nxticks,j=1:nxticks]
    for k=1:n 
        kmin,kmax = minimum(karr[k,:,:]),maximum(karr[k,:,:])
        ax = CairoMakie.Axis3(fig[k,1],
            xlabel = L"$x_1$",
            ylabel = L"$x_2$",
            zlabel = L"$K^{(\beta_1,\beta_2)}(x_1,\; x_2 \; |  \; \alpha, \; \gamma, \; \eta)",
            title = latexstring("\$\\beta_1 = $(β1[k]), \\; \\beta_2 = $(β2[k]), \\; \\alpha = $(α[k]), \\; \\gamma = $(γ[k]), \\; \\eta = $(η[k])\$"),
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
plot_gp_kernel_latticeseqb2_1s_contsurfs(;β1::Vector{Int64}=[0,1,0,1],β2::Vector{Int64}=[0,0,1,1],α::Vector{Int64}=[2,2,2,2],γ::Vector{Float64}=[1.,1.,1.,1.],η::Vector{Float64}=[1.,1.,1.,1.],xmin::Float64=-.1,xmax::Float64=1.1,nxticks::Int64=128,backgroundcolor::Symbol=:white) = plot_gp_kernel_1s_contsurfs(kernel_shiftinvar,β1,β2,α,γ,η,xmin,xmax,nxticks,backgroundcolor)
plot_gp_kernel_digitalseqb2g_1s_contsurfs(;β1::Vector{Int64}=[0,1,1],β2::Vector{Int64}=[0,0,1],α::Vector{Int64}=[4,4,4],γ::Vector{Float64}=[1.,1.,1.],η::Vector{Float64}=[1.,1.,1.],xmin::Float64=0.,xmax::Float64=1.,nxticks::Int64=128,backgroundcolor::Symbol=:white) = plot_gp_kernel_1s_contsurfs(kernel_digshiftinvar,β1,β2,α,γ,η,xmin,xmax,nxticks,backgroundcolor)

function plot_gp_optimization(gp::FastGaussianProcess;backgroundcolor::Symbol=:white)
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
    CairoMakie.xlims!(axζ,0,noptsp1-1)
    CairoMakie.linkxaxes!(axloss,axη)
    CairoMakie.linkxaxes!(axγ,axζ)
    return fig
end

function plot_gp_1s(gp::FastGaussianProcess;f::Union{Nothing,Function}=nothing,β::Vector{Int64}=[0],uncertainty::Float64=.05,xmin::Float64=0.,xmax::Float64=1.,nxticks::Int64=1024,markersize::Float64=16.,backgroundcolor::Symbol=:white)
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
        yhatticks = map(xtick->gp([xtick],[po]),xticks)
        CairoMakie.scatter!(ax,xticks,yhatticks,color=JULIA4LOGOCOLORS[1],markersize=markersize/2,label=latexstring("\$m_n^{($po)}(x)\$"))
        stdhatticks = sqrt.(map(xtick->var_post(gp,[xtick],[po]),xticks))
        ci_low,ci_high = yhatticks.-q*stdhatticks,yhatticks.+q*stdhatticks
        CairoMakie.band!(ax,xticks,ci_low,ci_high,color=(JULIA4LOGOCOLORS[1],.25),label=latexstring("\$m_n^{($po)}(x) \\pm $(round(q,digits=2)) \\; \\sigma_n^{($pi)}(x)\$"))
        if idx!==nothing CairoMakie.scatter!(ax,gp.x[:,1],gp.y[:,idx],markersize=markersize,color=:black,label=latexstring("\$(y^{($po)}_i)_{i=1}^{$(gp.n)}\$")) end 
        CairoMakie.Legend(fig[2*i-1,1],ax,orientation=:horizontal,framevisible=false) 
    end 
    fig 
end

function plot_gp_2s(gp::FastGaussianProcess;f::Union{Nothing,Function}=nothing,β::Matrix{Int64}=[0 0;],xmin::Float64=0.,xmax::Float64=1.,nxticks::Int64=32,markersize::Float64=16.,backgroundcolor::Symbol=:white)
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
        yhatticks = [gp([xticks[i],xticks[j]],po) for i=1:nxticks,j=1:nxticks]
        ax = CairoMakie.Axis3(fig[i,f===nothing ? 1 : 3],xlabel=L"$x_1$",ylabel=L"$x_2$",zlabel="",title=latexstring("\$m_n^{($po1,$po2)}(x)\$")); CairoMakie.xlims!(ax,xmin,xmax); CairoMakie.ylims!(ax,xmin,xmax)
        CairoMakie.surface!(ax,xticks,xticks,yhatticks,colormap=:julia_colorscheme)
        if idx!==nothing CairoMakie.scatter!(ax,gp.x[:,1],gp.x[:,2],gp.y[:,idx],markersize=markersize,color=:black,label=latexstring("\$(y^{($po1,$po2)}_i)_{i=1}^{$(gp.n)}\$")) end 
        ax = CairoMakie.Axis(fig[i,f===nothing ? 2 : 4],xlabel=L"$x_1$",ylabel=L"$x_2$",aspect=1,title=latexstring("\$m_n^{($po1,$po2)}(x)\$")); CairoMakie.xlims!(ax,xmin,xmax); CairoMakie.ylims!(ax,xmin,xmax)
        CairoMakie.heatmap!(ax,xticks,xticks,yhatticks,colormap=:julia_colorscheme)
        if idx!==nothing CairoMakie.scatter!(ax,gp.x[:,1],gp.x[:,2],markersize=markersize,color=:black,label=latexstring("\$(y^{($po1,$po2)}_i)_{i=1}^{$(gp.n)}\$")) end 
    end
    fig
end 
