module Nuclei
using Printf
using ..Basic

export Nucleus, RPAinfo, Set_nucleus, State
export Write_levels, Write_density
#
#
#
mutable struct Nucleus
    A::Int
    Z::Int
    N::Int
    R::Float64
    nmax
    eh
    jh
    lh
    node
    psi
    delta
    eqp
    v2
    eferm
    radii
    nuclmassc
    rho
    tau
    sj
    rhopair
    Etot::Float64
    Epair
    Binv0::Float64
    dpsi
    ddpsi
    Nucleus()=new(0, 0, 0, 0.0, [0,0],
                  zeros(Float64,nstatemax,2),
                  zeros(Int64,nstatemax,2), zeros(Int64,nstatemax,2),
                  zeros(Int64,nstatemax,2),
                  zeros(Float64,ngrid,nstatemax,2),
                  zeros(Float64,nstatemax,2), zeros(Float64,nstatemax,2),
                  zeros(Float64,nstatemax,2),
                  [-8.0, -8.0],
                  zeros(Float64,3),
                  nuclmass,
                  zeros(Float64,ngrid,2), zeros(Float64,ngrid,2),
                  zeros(Float64,ngrid,2), zeros(Float64,ngrid,2),
                  0.0, [0.0,0.0], 0.0, 
                  zeros(Float64,ngrid,nstatemax,2),
                  zeros(Float64,ngrid,nstatemax,2))
end
#
#
#
function Set_nucleus(Z::Int64, A::Int64)
    nucl=Nucleus()
    nucl.Z=Z
    nucl.A=A
    nucl.N=A-Z
    nucl.nuclmassc=nuclmass/(1.0-1.0/A)
    nucl.Binv0=hb2/2/nucl.nuclmassc
    nucl.R=1.2*(A-1)^(1.0/3.0)
    return nucl
end
#
#
#
function State(nucl::Nucleus, i::Int64, iq::Int64)
    return nucl.eh[i,iq], nucl.jh[i,iq], nucl.lh[i,iq],nucl.psi[:,i,iq],sqrt(nucl.v2[i,iq]),sqrt(1.0-nucl.v2[i,iq])
end
#
#
#
function Write_psi(nucl::Nucleus)
    f=open("single_wavefunction.out","w")
    for iq=1:2
        for i=1:nucl.nmax[iq]
            for ir=1:ngrid
                @printf(f,"%7.2f %12.5e  %2d %2d %9.4f %2d %2d\n",
                        ir*dr, nucl.psi[ir,i,iq], i, iq,
                        nucl.eh[i,iq],nucl.jh[i,iq],nucl.lh[i,iq])
            end
            @printf(f,"\n")
        end
    end
    close(f)
end
#
#
#
function Write_levels(nucl::Nucleus)
    f=open("single_particle.out","w")
    println(  "neutron")
    println(  " Num      energy   node j  l   occupation")
    println(f,"neutron")
    println(f, " Num      energy   node j  l   occupation")
    for i=nucl.nmax[2]:-1:1
        @printf("%4d  %12.7f  %2d %2d %2d %12.7f\n",
                i, nucl.eh[i,2], nucl.node[i,2],
                nucl.jh[i,2],nucl.lh[i,2],nucl.v2[i,2])
        @printf(f,"%4d  %12.7f  %2d %2d %2d %12.7f\n",
                i, nucl.eh[i,2], nucl.node[i,2],
                nucl.jh[i,2],nucl.lh[i,2],nucl.v2[i,2])
    end
    @printf(  "FermiE%12.7f\n",nucl.eferm[2])
    @printf(f,"FermiE%12.7f",nucl.eferm[2])
    #
    println(  " ")
    println(  "proton")
    println(   " Num      energy   node j  l   occupation")
    println(f," ")
    println(f,"proton")
    println(f, " Num      energy   node j  l   occupation")
    for i=nucl.nmax[1]:-1:1
        @printf("%4d  %12.7f  %2d %2d %2d %12.7f\n",
                i, nucl.eh[i,1], nucl.node[i,1],
                nucl.jh[i,1],nucl.lh[i,1],nucl.v2[i,1])
        @printf(f,"%4d  %12.7f  %2d %2d %2d %12.7f\n",
                i, nucl.eh[i,1], nucl.node[i,1],
                nucl.jh[i,1],nucl.lh[i,1],nucl.v2[i,1])
    end
    @printf(  "FermiE%12.7f\n",nucl.eferm[1])
    @printf(f,"FermiE%12.7f",nucl.eferm[1])
    close(f)
end
#
#
#
function Write_density(nucl::Nucleus)
    f=open("density.out","w")
    println(f,"      r   proton rho   neutron rho")
    for ir=1:ngrid
        r=ir*dr
        @printf(f,"%6.3f %12.7f %12.7f\n",r,nucl.rho[ir,1],nucl.rho[ir,2])
    end
    close(f)
end

end
