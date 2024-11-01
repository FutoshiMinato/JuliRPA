include("basic.jl")
include("math.jl")
include("nuclei.jl")
include("shfbcs.jl")
include("crpa.jl")
#include("angle.jl")

#using SymPy
using .Basic
using .Nuclei
using .Shfbcs
using .Math
using .Crpa
#using .Angle

function main()
    t0,t1,t2,t3,x0,x1,x2,x3,w0,w0p,alp,vpair,ij2=SkyrmeParameter("sko")
    t0,t1,t2,t3,x0,x1,x2,x3,w0,w0p,alp,vpair,ij2=SkyrmeParameter("sno")
    skyrme=Skyrme(t0,t1,t2,t3,x0,x1,x2,x3,w0,w0p,alp,vpair,ij2)
    #Nucleus
    O16 =Set_nucleus(8, 16)
    #SHFBCS
    mf=MeanField()
    fill!(mf.Binv, O16.Binv0)
    Calc_shfbcs(O16, skyrme, mf)
    Calc_continuum(O16, skyrme, mf)
end

@time main()

