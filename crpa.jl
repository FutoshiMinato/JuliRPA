module Crpa
using Printf
using CGcoefficient
using ..Basic
using ..Nuclei
using ..Shfbcs
using ..Math

export Calc_continuum, CRPA

mutable struct CRPA
    de
    ex
    str_free
    str_rpa
    L
    ISV
    scalingf
    gam
    CRPA()=new(0.2, 
        Vector{Float64}(undef, 0),
        Vector{Float64}(undef, 0),
        Vector{Float64}(undef, 0),
        1,
        0,
        1.0,
        1.0
        )
end

function Calc_continuum_Ex(nucl::Nucleus, skm::Skyrme, mf::MeanField, ex, gam, ISV, L, scalingf)
    LL=2*L
    vres=vresidual(nucl, skm, scalingf)
    vext=externalfield(L,ISV)
    G=zeros(Complex{Float64},2*ngrid,2*ngrid)
    A=zeros(Complex{Float64},2*ngrid,2*ngrid)
    for iq=1:2
        if iq == 1
            mr=0
        else
            mr=ngrid
        end
        for i=1:nucl.nmax[iq]
            if nucl.eh[i,iq] > nucl.eferm[iq]
    	        continue
            end
	    psi=nucl.psi[:,i,iq]
	    eh=nucl.eh[i,iq]
	    jh=nucl.jh[i,iq]
	    lh=nucl.lh[i,iq]
            Jmin=abs(LL-jh)
    	    Jmax=    LL+jh
	    for Jp=Jmin:Jmax
	        if mod((Jp+1)/2+lh+L,2) != 0
	            Lp=(Jp-1)/2
	        else
	            Lp=(Jp+1)/2
	        end
	        temp=(Jp+1)*(jh+1)/(2*L+1)*float(CG(jh//2, Jp//2, L, 1//2, -1//2, 0))^2*pi4m
	        up ,wp =green(mf, nucl,  ex+eh,  gam, iq, Lp, Jp)
	        up1,wp1=green(mf, nucl, -ex+eh, -gam, iq, Lp, Jp)
	        for ir=1:ngrid
		    for jr=ir:ngrid
	  	        G[ir+mr,jr+mr]+=(temp*psi[ir]*psi[jr]
			*(up[ir]*wp[jr]+up1[ir]*wp1[jr]))
		    end
		end
            end
	end
    end
    #
    for ir=1:2*ngrid
        for jr=1:2*ngrid
            G[jr,ir]=G[ir,jr]
            if ir <= ngrid
	        if jr <= ngrid
	            A[ir,jr]=G[ir,jr      ]*vres[jr      ,2]
		    #G11*VIS
	        else
	            A[ir,jr]=G[ir,jr-ngrid]*vres[jr-ngrid,1]
		    #G11*VIV
	       end
	    else
	        if jr <= ngrid
	            A[ir,jr]=G[ir,jr+ngrid]*vres[jr      ,1]
		    #G22*VIV
	        else
	            A[ir,jr]=G[ir,jr      ]*vres[jr-ngrid,2]
		    #G22*VIS
	        end
	    end
	    if ir == jr
	        A[ir,jr]+=1.0
	    end
        end
    end
    #
    B=inv(A)
    #
    drhof,drhorpa,sfree,srpa=Calc_Strength(G,B,vext)
    return drhof,drhorpa,sfree,srpa
end
#
#
#
function Calc_Strength(G,B,vext)
    #
    sfree=0.0
    drhof=zeros(Complex{Float64},2*ngrid)
    for ir=1:2*ngrid
    	for jr=1:2*ngrid
	    drhof[ir]+=G[ir,jr]*vext[jr]*dr
	end
	sfree+=vext[ir]*drhof[ir]*dr
    end
    #
    srpa=0.0
    drhorpa=zeros(Complex{Float64},2*ngrid)
    for ir=1:2*ngrid
    	for jr=1:2*ngrid
	    drhorpa[ir]+=B[ir,jr]*drhof[jr]
	end
	srpa+=vext[ir]*drhorpa[ir]*dr
    end
    return drhof,drhorpa,sfree,srpa
end
#
#
#
function vresidual(nucl::Nucleus, skm::Skyrme, scalingf)
    t0=skm.t0
    x0=skm.x0
    t3=skm.t3
    x3=skm.x3
    alp=skm.alp
    vres=zeros(Float64,ngrid,2)
    for ir=1:ngrid
        r=ir*dr
    	rhot=nucl.rho[ir,2]+nucl.rho[ir,1]
    	rhou=nucl.rho[ir,2]-nucl.rho[ir,1]
    	rhoa=rhot^alp
    	rhoa1=rhot^(alp-1)
        v01= (0.75*t0 
             +0.0625*t3*(alp+2)*(alp+1)*rhoa 
             -0.0208333*t3*(1+2*x3)*alp*(alp-1)*rhou^2/rhot^2*rhoa)
        v02= (-0.25*t0*(1+2*x0)-0.0416666*t3*(1+2*x3)*rhoa)
        v01r=-t3/24*(1+2*x3)*alp*rhoa1*rhou
        v02r=-t3/24*(1+2*x3)*alp*rhoa1*rhou
    	vres[ir,1]=scalingf*(v01+v01r-v02-v02r)/r^2*dr # isovector V-V_tau
    	vres[ir,2]=scalingf*(v01+v01r+v02+v02r)/r^2*dr # isoscalar V+V_tau
    end
    return vres
end
#
#
#
function externalfield(L,ISV)
    f=zeros(Float64,2*ngrid)
    K=L
    if L == 0
        K=2
    end
    for ir=1:ngrid
        r=ir*dr
        f[ir]=r^K
        if ISV == 0
            f[ir+ngrid]= f[ir]
        elseif ISV == 1
            f[ir+ngrid]=-f[ir]
        end
    end
    return f
end

function Vcent(l, r, nucl::Nucleus)
    return l*(l+1)*hb2/2/nucl.nuclmassc/r^2
end

function green(mf::MeanField, nucl::Nucleus, e, gam, iq, l, j)
    up=zeros(Complex{Float64},ngrid)
    wp=zeros(Complex{Float64},ngrid)
    S =zeros(Complex{Float64},ngrid)
    V=zeros(Float64,ngrid)
    fac=dr^2*(2*nucl.nuclmassc/hb2)
    xls=(j/2.0)*(j/2.0+1.0)-l*(l+1)-3.0/4.0
    for ir=1:ngrid
        r=ir*dr
        V[ir]=(mf.Vn[ir,iq]+mf.Vm[ir,iq]+mf.Vls[ir,iq]*xls
        +(mf.Vc[ir]+mf.Vcex[ir])*(2-iq)
        +Vcent(l,r,nucl)*mf.Binv[ir,iq]/nucl.Binv0)
    end
  
    ecomp=e+gam*im/2

    psi0=dr^(l+1)
    psi1=(2+5.0/6.0*fac*nucl.Binv0/mf.Binv[1,iq]*(V[1]-ecomp))*psi0
    psi1=psi1/(1-fac/12*nucl.Binv0/mf.Binv[2,iq]*(V[2]-ecomp))

    up[1]=psi0
    up[2]=psi1

    for ir=3:ngrid
        dd=psi0-2*psi1
        dd-=fac/12*nucl.Binv0/mf.Binv[ir-2,iq]*(V[ir-2]-ecomp)*psi0
        dd-=5*fac/6*nucl.Binv0/mf.Binv[ir-1,iq]*(V[ir-1]-ecomp)*psi1
        S[ir]=fac*nucl.Binv0/mf.Binv[ir-2,iq]*(V[ir-2]-ecomp)

        cc=1-fac/12*nucl.Binv0/mf.Binv[ir,iq]*(V[ir]-ecomp)
        cin=1.0/cc

        up[ir]=-cin*dd

        psi0=psi1
        psi1=up[ir]
    end
  
  #initialize irregular solution to outgoing wave
    ZZ=-real(S[ngrid-1])
    if ZZ > 0
        ak=sqrt(ZZ)
        wp[ngrid-1]=1.0
        wp[ngrid]=(1.0-ak^2/2.0)+im*ak
    else
        wp[ngrid-1]=0.000
        wp[ngrid]  =0.001
    end
    psi0=wp[ngrid]
    psi1=wp[ngrid-1]

    for ir=ngrid-2:-1:1
        dd=psi0-2*psi1
        dd-=fac/12*nucl.Binv0/mf.Binv[ir+2,iq]*(V[ir+2]-ecomp)*psi0
        dd-=5*fac/6*nucl.Binv0/mf.Binv[ir+1,iq]*(V[ir+1]-ecomp)*psi1

        cc=1-fac/12*nucl.Binv0/mf.Binv[ir,iq]*(V[ir]-ecomp)
        cin=1.0/cc

        wp[ir]=-cin*dd

        psi0=psi1
        psi1=wp[ir]
    end

    for ir=1:ngrid
        red=sqrt(nucl.Binv0/mf.Binv[ir,iq])
        cc=up[ir]*red
        dd=wp[ir]*red
        up[ir]=cc
        wp[ir]=dd
    end

    Rn=1.5*nucl.A^(1.0/3.0)

    is=Int64(floor(Rn/dr))
    dup0=up[is-2]-8.0*up[is-1]+8.0*up[is+1]-up[is+2]
    dup0=dup0/(12*dr)
    dwp0=wp[is-2]-8.0*wp[is-1]+8.0*wp[is+1]-wp[is+2]
    dwp0=dwp0/(12*dr)
    ron=(wp[is]*dup0-up[is]*dwp0)*mf.Binv[is,iq]

    for ir=1:ngrid
        wp[ir]=wp[ir]/ron
    end

  return up,wp
end
#
#
#
function Calc_continuum(nucl::Nucleus, skm::Skyrme, mf::MeanField, crpa::CRPA)
    de=crpa.de
    f=open("strength.out","w")
    sftot=0.0
    srtot=0.0
    for i=0:200
        ex=i*de
        drhof,drhorpa,sfree,srpa=Calc_continuum_Ex(nucl, skm, mf, ex, crpa.gam, crpa.ISV, crpa.L, crpa.scalingf)
#        @printf(  "%8.3f %12.7f %12.7f\n", ex, imag(sfree)/pi, imag(srpa)/pi)
        @printf(f,"%8.3f %12.7f %12.7f\n", ex, imag(sfree)/pi, imag(srpa)/pi)
        append!(crpa.ex,ex)
        append!(crpa.str_free,imag(sfree)/pi)
        append!(crpa.str_rpa ,imag(srpa)/pi)
    	sftot+=ex*imag(sfree)/pi*de
    	srtot+=ex*imag(srpa) /pi*de
    end
    close(f)
#    @printf("         %12.7f %12.7f\n", sftot, srtot)
end

end