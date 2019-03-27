using ESDLPlots
using Test

# write your own tests here
using ESDL
using Widgets
using Compose
import Plots

c = Cube()
d = subsetcube(c,variable = ["gross_primary_productivity","net_ecosystem_exchange","terrestrial_ecosystem_respiration]"], time=Date(2001)..Date(2001,1,15))
pm = plotMAP(d)
pm1 = plotMAP(d,variable="gross",time=Date(2001))
@test isa(pm,Widget)
@test isa(pm1,Context)

pts = plotTS(d);

pts1 = plotTS(d,lon=35,lat=55,var="net");


@test isa(pts,Widget)
@test isa(pts1,Plots.Plot)

pxy = plotXY(d,xaxis="var")

pscatter = plotScatter(d)

pimage = plotMAPRGB(d)

@test isa(pimage,Widget)
@test isa(pxy, Widget)
@test isa(pscatter, Widget)
