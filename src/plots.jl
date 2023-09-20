function plot_gp_optimization(gp::GaussianProcessLatticeSeqB2)
    noptsp1 = length(gp.losses)
    @assert noptsp1>1
    xrange = [i for i=0:noptsp1-1]
    fig = CairoMakie.Figure(resolution=(1000,700),sharex=true)
    axloss = CairoMakie.Axis(fig[1,1],
        ylabel = L"$\text{Loss } L(\theta | Y)$")
    CairoMakie.lines!(axloss,xrange,gp.losses,color=JULIA4LOGOCOLORS[1],linewidth=3)
    CairoMakie.xlims!(axloss,0,noptsp1-1)
    axγ = CairoMakie.Axis(fig[1,2],
        ylabel = L"$\gamma$",
        yscale = log10)
    CairoMakie.lines!(axγ,xrange,gp.γs,color=JULIA4LOGOCOLORS[2],linewidth=3)
    CairoMakie.xlims!(axγ,0,noptsp1-1)
    axη = CairoMakie.Axis(fig[2,1],
        xlabel = "step",
        ylabel = L"$\eta$",
        yscale = log10)
    for j=1:gp.s CairoMakie.lines!(axη,xrange,gp.ηs[:,j],color=(JULIA4LOGOCOLORS[3],(gp.s-j+1)/gp.s),linewidth=3) end 
    CairoMakie.xlims!(axη,0,noptsp1-1)
    axζ = CairoMakie.Axis(fig[2,2],
        xlabel = "step",
        ylabel = L"$\zeta$",
        yscale = log10)
    for l=1:gp.r CairoMakie.lines!(axζ,xrange,gp.ζs[:,l],color=(JULIA4LOGOCOLORS[4],(gp.r-l+1)/gp.r),linewidth=3) end 
    CairoMakie.xlims!(axζ,0,noptsp1-1)
    CairoMakie.linkxaxes!(axloss,axη)
    CairoMakie.linkxaxes!(axγ,axζ)
    return fig
end