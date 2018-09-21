mutable struct XYPlot <: ESDLPlot
  xaxis
  group
end
plotAxVars(p::XYPlot)=[FixedAx(p.xaxis,"X Axis",true,false,1),FixedAx(p.group,"Group",false,false,2)]

function plotCall(::XYPlot,d::AbstractCubeData,ixaxis,igroup,otherinds...)

  axlist = caxes(d)
  inds1 = ntuple(  i->in(i,(ixaxis,igroup)) ? (:) : axVal2Index(axlist[i],otherinds[i]), length(otherinds))
  x1 = d[inds1...]

  if igroup > 0
    igroup < ixaxis && (x1=transpose(x1))
    if isa(axlist[ixaxis],CategoricalAxis)
      plotf = StatPlots.groupedbar
      x1 = x1'
    else
      plotf = Plots.plot
    end
    labs = reshape(string.(axlist[igroup].values),(1,length(axlist[igroup])))
    xlabel = axname(axlist[ixaxis])
    p=plotf(axlist[ixaxis].values,x1,
    lab=labs,
    xlabel=xlabel)
  else
    plotf = isa(axlist[ixaxis],CategoricalAxis) ? Plots.bar : Plots.plot
    p=plotf(axlist[ixaxis].values,x1,xlabel=axname(axlist[ixaxis]),fmt=:png)
  end
  p
end
nplotCubes(::XYPlot)=1

"""
    plotXY(cube::AbstractCubeData; group=0, xaxis=-1, kwargs...)

Generic plotting tool for cube objects, can be called on any type of cube data.

### Keyword arguments

* `xaxis` which axis is to be used as x axis. Can be either an axis Datatype or a string. Short versions of axes names are possible as long as the axis can be uniquely determined.
* `group` it is possible to group the plot by a categorical axis. Can be either an axis data type or a string.
* `dim=value` can set other dimensions to certain values, for example `lon=51.5` will fix the longitude for the resulting plot

If a dimension is not the x axis or group variable and is not fixed through an additional keyword, a slider or dropdown menu will appear to select the axis value.
"""
function plotXY(cube::AbstractCubeData{T};group=nothing,xaxis=nothing,kwargs...) where T

  return plotGeneric(XYPlot(xaxis,group),cube;kwargs...)

end

mutable struct ScatterPlot <: ESDLPlot
  vsaxis
  alongaxis
  group
  c_1
  c_2
end
plotAxVars(p::ScatterPlot)=[
  FixedAx(p.vsaxis,"VS Axis",true,true,-1),
  FixedAx(p.alongaxis,"Along",true,false,1),
  FixedAx(p.group,"Group",false,false,2),
  FixedVar(p.vsaxis,p.c_1,"X Axis",true),
  FixedVar(p.vsaxis,p.c_2,"Y Axis",true)
  ]

function plotCall(::ScatterPlot,d::AbstractCubeData,ivsaxis,ialongaxis,igroup,c1,c2,otherinds...)

  axlist = caxes(d)
  inds1 = ntuple(  i->(i==ivsaxis)    ? axVal2Index(axlist[i],c1)  : (i==ialongaxis) ? (:) : (i==igroup)     ? (:) : axVal2Index(axlist[i],otherinds[i]), length(otherinds))
  inds2 = ntuple(  i->(i==ivsaxis)    ? axVal2Index(axlist[i],c2)  : (i==ialongaxis) ? (:) : (i==igroup)     ? (:) : axVal2Index(axlist[i],otherinds[i]), length(otherinds))
  x1 = d[inds1...]
  x2 = d[inds2...]

  goodinds = map((m1,m2)->!ismissing(m1) && !ismissing(m2),x1,x2)
  a_1 = x1[goodinds]
  a_2 = x2[goodinds]
  if isempty(a_1)
    push!(a_1,1);push!(a_2,1)
  end
  nPoints = length(a_1)
  pointSize = min(5000/nPoints,3)
  msw=pointSize > 2 ? 1 : 0
  fmt=nPoints>20000 ? :png : :svg
  if igroup > 0 && (igroup != ialongaxis)
    plotf = isa(axlist[ivsaxis],CategoricalAxis) ? StatPlots.groupedbar : Plots.plot
    igroup < ialongaxis && (a_1=transpose(a_1);a_2=transpose(a_2))
    p=plotf(a_1,a_2,
      xlabel=string(axlist[ivsaxis].values[inds1[ivsaxis]]),
      ylabel=string(axlist[ivsaxis].values[inds2[ivsaxis]]),
      lab=reshape(string.(axlist[igroup].values),(1,length(axlist[igroup]))),
      fmt=fmt,
      ms=pointSize,
      markerstrokewidth=msw,
      )
  else
    p=Plots.scatter(a_1,a_2,
      xlabel=string(axlist[ivsaxis].values[inds1[ivsaxis]]),
      ylabel=string(axlist[ivsaxis].values[inds2[ivsaxis]]),
      lab="",
      fmt=fmt,
      ms=pointSize,
      markerstrokewidth=msw,
      )
  end
  p
end
"""
    plotScatter(cube::AbstractCubeData; vsaxis=VariableAxis, alongaxis=0, group=0, xaxis=0, yaxis=0, kwargs...)

Generic plotting tool for cube objects to generate scatter plots, like variable `A` against variable `B`. Can be called on any type of cube data.

### Keyword arguments

* `vsaxis` determines the axis from which the x and y variables are drawn.
* `alongaxis` determines the axis along which the variables are plotted. E.g. if you choose `TimeAxis`, a dot will be plotted for each time step.
* `xaxis` index or value of the variable to plot on the x axis
* `yaxis` index or value of the variable to plot on the y axis
* `group` it is possible to group the plot by an axis. Can be either an axis data type or a string. *Caution: This will increase the number of plotted data points*
* `dim=value` can set other dimensions to certain values, for example `lon=51.5` will fix the longitude for the resulting plot

If a dimension is not the `vsaxis` or `alongaxis` or `group` and is not fixed through an additional keyword, a slider or dropdown menu will appear to select the axis value.
"""
function plotScatter(cube::AbstractCubeData{T};group=nothing,vsaxis=VariableAxis,xaxis=nothing,yaxis=nothing,alongaxis=nothing,kwargs...) where T

  return plotGeneric(ScatterPlot(vsaxis,alongaxis,group,xaxis,yaxis),cube;kwargs...)

end

export plotHist
import ESDL.Proc.DATOnlineStats: HistogramCube, tohist
plotHist(c::HistogramCube)=plotScatter(tohist(c),vsaxis="Hist",xaxis="MidPoints",yaxis="Frequency",alongaxis="Bin")
