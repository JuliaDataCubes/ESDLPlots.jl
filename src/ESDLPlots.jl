module ESDLPlots
export plotTS, plotMAP, plotXY, plotScatter, plotMAPRGB
export plotlyjs, gadfly, gr, pyplot
using ESDL.Cubes
using ESDL.CubeAPI
using ESDL.CubeAPI.Mask
import ESDL.DAT
import ESDL.DAT: findAxis,getFrontPerm
import ESDL.Cubes.Axes.axname
import Reactive: Signal
import Interact: slider, dropdown, Observable, observe
import Colors: RGB, @colorant_str, colormap,  distinguishable_colors
import FixedPointNumbers: Normed
import Base.Cartesian: @ntuple,@nexprs
import Measures
import Compose
import Images
import DataStructures: OrderedDict
import StatPlots
import PlotUtils: optimize_ticks, cgrad
import Compose: rectangle, text, line, compose, context, stroke, svgattribute, bitmap, HCenter, VBottom, HRight, VCenter

function __init()__
  eval(:(import Plots))
  eval(:(import Plots: plotlyjs, gr, pyplot))
end

const U8=Normed{UInt8,8}


abstract type ESDLPlot end
"Expression to evaluate after the data is loaded"
getafterEx(::ESDLPlot)=Expr(:block)

"Setting fixed variables"
getFixedVars(::ESDLPlot,cube)=Expr(:block)

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


toYr(tx::TimeAxis)=((tx.values.startyear+(tx.values.startst-1)/tx.values.NPY):(1.0/tx.values.NPY):(tx.values.stopyear+(tx.values.stopst-1)/tx.values.NPY))-(tx.values.startyear+(tx.values.startst-1)/tx.values.NPY)

r1(x)=reshape(x,length(x))
prepAx(x)=x.values
prepAx(x::TimeAxis)=toYr(x)
function repAx(x,idim,ax)
  l=length(x)
  inrep=prod(size(x)[1:idim-1])
  outrep=div(l,(inrep*size(x,idim)))
  repeat(collect(ax),inner=[inrep],outer=[outrep])
end
function count_to(f,c,i)
  ni=0
  for ind=1:i
    f(c[ind]) && (ni+=1)
  end
  return ni
end

getWidget(x::CategoricalAxis;label=axname(x))       = dropdown(Dict(zip(x.values,1:length(x.values))),label=label)
getWidget(x::RangeAxis{T};label=axname(x)) where {T<:Real} = step(x.values) > 0 ? slider(x.values,label=label) : slider(reverse(x.values),label=label)
getWidget(x::RangeAxis;label=axname(x))             = slider(x.values,label=label)
getWidget(x::SpatialPointAxis;label="Spatial Point")= slider(1:length(x),label=label)

plotTS(x;kwargs...)=plotXY(x;xaxis=TimeAxis,kwargs...)

function setPlotAxis(a::FixedAx,axlist,fixedAxes,customobs,positionobs)
  ix=a.axis==nothing ? 0 : findAxis(a.axis,axlist)
  if ix>0
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
    push!(customobs,axVal2Index(axlist[a.depAxis],a.varVal))
    positionobs[a.depAxis]=0
  else
    push!(customobs,a)
  end
end

import Interact: observe, Widget

function createWidgets(axlist,availableAxis,availableIndices,axlabels,widgets,axtuples,customobs,positionobs)

  if !isempty(availableAxis)
    for (icust,at) in enumerate(customobs)
      if isa(at,FixedAx)
        options = collect(at.musthave ? zip(axlabels[availableIndices],availableIndices) : zip(["None";axlabels[availableIndices]],[0;availableIndices]))
        axmenu  = dropdown(OrderedDict(options),label=at.widgetlabel,value=options[1][2],value_label=options[1][1])
        sax = observe(axmenu)
        widgets[Symbol(at.widgetlabel)]=axmenu
        customobs[icust] = sax
      elseif isa(at,FixedVar)
        w=getWidget(axlist[at.depAxis],label=at.widgetlabel)
        widgets[Symbol(at.widgetlabel)]=w
        customobs[icust] = observe(w)
      else
        #println("Skipping selected")
      end
    end
    for i in availableIndices
      w=getWidget(axlist[i])
      widgets[Symbol(axname(axlist[i]))]=w
      positionobs[i] = observe(w)
    end
  else
    for (icust,at) in enumerate(customobs)
      if isa(at,FixedAx)
        at.musthave && error("No axis left to put on $label")
      end
    end
  end
end

const namedcolms=Dict(
:viridis=>[cgrad(:viridis)[ix] for ix in linspace(0,1,100)],
:magma=>[cgrad(:magma)[ix] for ix in linspace(0,1,100)],
:inferno=>[cgrad(:inferno)[ix] for ix in linspace(0,1,100)],
:plasma=>[cgrad(:plasma)[ix] for ix in linspace(0,1,100)])
typed_dminmax(::Type{T},dmin,dmax) where {T<:Integer}=(Int(dmin),Int(dmax))
typed_dminmax(::Type{T},dmin,dmax) where {T<:AbstractFloat}=(Float64(dmin),Float64(dmax))
typed_dminmax2(::Type{T},dmin,dmax) where {T<:Integer}=(isa(dmin,Tuple) ? (Int(dmin[1]),Int(dmin[2]),Int(dmin[3])) : (Int(dmin),Int(dmin),Int(dmin)), isa(dmax,Tuple) ? (Int(dmax[1]),Int(dmax[2]),Int(dmax[3])) : (Int(dmax),Int(dmax),Int(dmax)))
typed_dminmax2(::Type{T},dmin,dmax) where {T<:AbstractFloat}=(isa(dmin,Tuple) ? (Float64(dmin[1]),Float64(dmin[2]),Float64(dmin[3])) : (Float64(dmin),Float64(dmin),Float64(dmin)), isa(dmax,Tuple) ? (Float64(dmax[1]),Float64(dmax[2]),Float64(dmax[3])) : (Float64(dmax),Float64(dmax),Float64(dmax)))

mygetval(i)=i
mygetval(i::Observable)=i[]


import Interact
import InteractBase
import Observables
function plotGeneric(plotObj::ESDLPlot, cube::CubeAPI.AbstractCubeData{T};kwargs...) where T


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
      if ix > 0
        push!(fixedAxes,axlist[ix])
        positionobs[ix] = axVal2Index(axlist[ix],val,fuzzy=true)
      end
    end

  availableIndices=findall(ax->!in(ax,fixedAxes),axlist)
  availableAxis=axlist[availableIndices]

  createWidgets(axlist,availableAxis,availableIndices,axlabels,widgets,pAxVars,customobs,positionobs)

  ofirst = plotCall(plotObj,cube,map(mygetval,customobs)...,map(mygetval,positionobs)...);
  s=Observable(ofirst)

  map!((i...)->plotCall(plotObj,cube,i...),s,customobs...,positionobs...);

  #map(display,widgets)
  #display(s)
  layout = t -> InteractBase.node(:div, map(InteractBase.center, InteractBase.values(InteractBase.components(t)))..., map(InteractBase.center, s))
  Widget{:ESDLPlot}(widgets,output=Observables.throttle(1.0,s),layout = layout)
end

end # module
