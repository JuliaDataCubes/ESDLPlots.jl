module ESDLPlots
using ImageMagick
export plotTS, plotMAP, plotXY, plotScatter, plotMAPRGB
export plotlyjs, gadfly, gr, pyplot
using YAXArrays: Cubes
using YAXArrays.Cubes.Axes: axname, findAxis, CategoricalAxis, RangeAxis, CubeAxis,
  axVal2Index
using YAXArrays.Cubes: caxes
using Widgets: slider, dropdown, Widget, Widgets
using Observables: Observable, observe, throttle
using Colors: RGB, @colorant_str, colormap,  distinguishable_colors
using FixedPointNumbers: Normed
using Measures: Measures
using Compose: Compose, rectangle, text, line, compose, context,
  stroke, svgattribute, bitmap, HCenter, VBottom, HRight, VCenter
using Images: Images
using DataStructures: OrderedDict
using PlotUtils: optimize_ticks, cgrad
using Plots: Plots, plotlyjs, gr, pyplot, plot, bar, scatter
using StatsPlots: groupedbar

const U8=Normed{UInt8,8}


abstract type ESDLPlot end

mutable struct FixedAx
  axis
  widgetlabel::String
  musthave::Bool
  isimmu::Bool
  position::Int
end

mutable struct FixedVar
  depAxis
  varVal
  widgetlabel::String
  musthave::Bool
end

include("maps.jl")
include("other.jl")

getWidget(x::CategoricalAxis;label=axname(x))       = dropdown(x.values,label=label)
getWidget(x::RangeAxis{T};label=axname(x)) where {T<:Real} = (last(x.values)-first(x.values)) > 0 ? slider(x.values,label=label) : slider(reverse(x.values),label=label)
getWidget(x::RangeAxis;label=axname(x))             = slider(x.values,label=label)

plotTS(x;kwargs...)=plotXY(x;xaxis="time",kwargs...)

function setPlotAxis(a::FixedAx,axlist,fixedAxes,customobs,positionobs)
  ix = a.axis===nothing ? nothing : findAxis(a.axis,axlist)
  if ix !== nothing
    push!(fixedAxes,axlist[ix])
    push!(customobs,ix)
    a.axis=ix
    positionobs[ix]=0
  else
    push!(customobs,a)
    a.axis = a.isimmu ? error("Axis $a must be selected.") : 0
  end
end
function setPlotAxis(a::FixedVar,axlist,fixedAxes,customobs,positionobs)
  a.depAxis=findAxis(a.depAxis,axlist)
  if a.varVal!=nothing
    push!(customobs,a.varVal)
    positionobs[a.depAxis]=0
  else
    push!(customobs,a)
  end
end

cart(i::Integer) = CartesianIndex((i,))

function createWidgets(axlist,availableAxis,availableIndices,axlabels,widgets,axtuples,customobs,positionobs)

  if !isempty(availableAxis)
    for (icust,at) in enumerate(customobs)
      if isa(at,FixedAx)
        options = collect(at.musthave ? zip(axlabels[availableIndices],availableIndices) : zip(["None",axlabels[availableIndices]...],[0,availableIndices...]))
        axmenu  = dropdown(OrderedDict(options),label=at.widgetlabel,value=options[1][2],value_label=options[1][1])
        sax = observe(axmenu)
        widgets[Symbol(at.widgetlabel)]=axmenu
        customobs[icust] = sax
      elseif isa(at,FixedVar)
        w=getWidget(axlist[at.depAxis],label=at.widgetlabel)
        widgets[Symbol(at.widgetlabel)]=w
        customobs[icust] = throttle(1.0,observe(w))
      else
        #println("Skipping selected")
      end
    end
    for i in availableIndices
      w=getWidget(axlist[i])
      widgets[Symbol(axname(axlist[i]))]=w
      positionobs[i] = throttle(1.0,observe(w))
    end
  else
    for (icust,at) in enumerate(customobs)
      if isa(at,FixedAx)
        at.musthave && error("No axis left to put on $(at.widgetlabel)")
        customobs[icust] = -1
      end
    end
  end
end

const namedcolms=Dict(
:viridis=>[cgrad(:viridis)[ix] for ix in range(0,stop=1,length=100)],
:magma=>[cgrad(:magma)[ix] for ix in range(0,stop=1,length=100)],
:inferno=>[cgrad(:inferno)[ix] for ix in range(0,stop=1,length=100)],
:plasma=>[cgrad(:plasma)[ix] for ix in range(0,stop=1,length=100)])
typed_dminmax(::Type{<:Union{T,Missing}},dmin,dmax) where T = typed_dminmax(T,dmin,dmax)
typed_dminmax2(::Type{<:Union{T,Missing}},dmin,dmax) where T = typed_dminmax2(T,dmin,dmax)
typed_dminmax(::Type{T},dmin,dmax) where {T<:Integer}=(Int(dmin),Int(dmax))
typed_dminmax(::Type{T},dmin,dmax) where {T<:AbstractFloat}=(Float64(dmin),Float64(dmax))
typed_dminmax2(::Type{T},dmin,dmax) where {T<:Integer}=(isa(dmin,Tuple) ? (Int(dmin[1]),Int(dmin[2]),Int(dmin[3])) : (Int(dmin),Int(dmin),Int(dmin)), isa(dmax,Tuple) ? (Int(dmax[1]),Int(dmax[2]),Int(dmax[3])) : (Int(dmax),Int(dmax),Int(dmax)))
typed_dminmax2(::Type{T},dmin,dmax) where {T<:AbstractFloat}=(isa(dmin,Tuple) ? (Float64(dmin[1]),Float64(dmin[2]),Float64(dmin[3])) : (Float64(dmin),Float64(dmin),Float64(dmin)), isa(dmax,Tuple) ? (Float64(dmax[1]),Float64(dmax[2]),Float64(dmax[3])) : (Float64(dmax),Float64(dmax),Float64(dmax)))

mygetval(i)=i
mygetval(i::Observable)=i[]
mygetval(i::FixedAx)=-1

function plotGeneric(plotObj::ESDLPlot, cube;kwargs...)


  axlist=caxes(cube)

  axlabels=map(axname,axlist)
  fixedAxes=CubeAxis[]
  customobs=[]
  widgets=OrderedDict{Symbol,Any}()
  positionobs=Array{Any}(undef,ndims(cube))

  pAxVars=plotAxVars(plotObj)

  foreach(t->setPlotAxis(t,axlist,fixedAxes,customobs,positionobs),pAxVars)

  for (sy,val) in kwargs
    ix = findAxis(string(sy),axlist)
    if ix !== nothing
      push!(fixedAxes,axlist[ix])
      positionobs[ix] = val
    end
  end

  availableIndices=findall(ax->!in(ax,fixedAxes),axlist)
  availableAxis=axlist[availableIndices]

  createWidgets(axlist,availableAxis,availableIndices,axlabels,widgets,pAxVars,customobs,positionobs)

  allobs = [customobs;positionobs]

  if any(i->isassigned(allobs,i) && isa(allobs[i],Observable),1:length(allobs))

    s=Observable{Any}(" ")

    map!((i...)->plotCall(plotObj,cube,i...),s,customobs...,positionobs...);

    ifirstwig = findfirst(i->isassigned(allobs,i) && isa(allobs[i],Observable),1:length(allobs))
    firstwig = allobs[ifirstwig]
    firstwig[] = firstwig[]

    layout = (Widgets.manipulatelayout)((Widgets.get_backend)())
    return Widget{:ESDLPlot}(widgets,output=s,layout = layout)

  else

    return plotCall(plotObj,cube,customobs...,positionobs...)

  end
end

end # module
