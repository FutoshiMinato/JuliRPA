#
# Skyrme-Hartree-Fock (Julia)
#
# History: 2023/10/06
#
module Shfbcs
using Printf
using ..Nuclei
using ..Basic
using ..Math

export SkyrmeParameter, Skyrme, Get_param_skyrme, Calc_shfbcs
export MeanField, Write_potential

#iprint=true
iprint=false
rgrid=zeros(Float64,ngrid)
for ir=1:ngrid
    rgrid[ir]=ir*dr
end

function SkyrmeParameter(name)
    if name == "skm"
        return -2645.000,  410.000,-135.0000, 15595.000,
	        0.090000, 0.000000, 0.000000, 0.0000000,
		130.0000, 130.0000, 1.0/6.0, [0,0],
		false
    elseif name == "sko"
        return -2103.653,  303.352, 791.6740, 13553.252,
	       -0.210701,-2.810752,-1.461595,-0.4298810,
		353.1560,-397.4980, 0.25,    [0,0],
		false
    elseif name == "sko2"
        return -2099.419,  301.531, 154.781, 13526.464,
	       -0.029503,-1.325732,-2.323439,-0.147404,
		287.79,-165.7776, 0.25,    [0,0],
		false
    elseif name == "sly4"
        return -2488.91,  486.82, -546.39, 13777.0,
	        0.834,-0.3438,-1.0,1.3539,
		123.0,123.0, 0.1666666,    [0,0],
		false
    elseif name == "sly5"
        return -2484.88,  483.13, -549.4, 13763.0,
	        0.778,-0.328,-1.0, 1.267,
		126.0,126.0, 0.1666666,    [0,0],
		false
    elseif name == "sg2"
        return -2645.0,  340.0, -41.9, 15595.0,
	        0.09,-0.588, 1.425, 0.06044,
		105.0,105.0, 0.1666666,    [0,0],
		false
    elseif name == "sno"
        return -1100.000,  0.0, 0.0, 15000.000,
	        0.0,      0.0,      0.0,      0.0,
		0.0,      0.0,      1.0,     [0,0],
		false
    else
        println("No Skyrme Force requested!")
    end
end

struct Skyrme
    t0::Float64
    t1::Float64
    t2::Float64
    t3::Float64
    x0::Float64
    x1::Float64
    x2::Float64
    x3::Float64
    w0::Float64
    w0p::Float64
    alp::Float64
    vpair::Array{Float64,1}
    ij2::Bool
end

function Get_param_skyrme(skm::Skyrme)
    return skm.t0,skm.t1,skm.t2,skm.t3,skm.x0,skm.x1,skm.x2,skm.x3,skm.w0,skm.w0p,skm.alp
end

mutable struct MeanField
    V
    Vn
    Vls
    Vm
    Vc
    Vcex
    DeltaPot
    Binv
    MeanField()=new(zeros(Float64, ngrid, 2), zeros(Float64, ngrid, 2),
                    zeros(Float64, ngrid, 2), zeros(Float64, ngrid, 2),
                    zeros(Float64, ngrid), zeros(Float64, ngrid),
                    zeros(Float64, ngrid, 2),
                    zeros(Float64, ngrid, 2))
end

function WoodSaxon(mf::MeanField, nucl::Nucleus)
    V0  =-58.0
    Vls0= 15.0
    adiff=0.65
    for ir=1:ngrid
        r=rgrid[ir]
        f =1.0/(1+exp((r-nucl.R)/adiff))
        df=-1.0/(1.0+exp((r-nucl.R)/adiff))^2*exp((r-nucl.R)/adiff)/adiff

        fac=1.0-0.67*(nucl.N-nucl.Z)/nucl.A
        mf.Vn[ir,2]  =fac*V0*f
        mf.Vls[ir,2] =fac*Vls0*df/r
        mf.Vm[ir,2]=0.0
        mf.DeltaPot[ir,2]=12.0/sqrt(nucl.A)

        fac=1.0+0.67*(nucl.N-nucl.Z)/nucl.A
        mf.Vn[ir,1]  =fac*V0*f
        mf.Vls[ir,1] =fac*Vls0*df/r
        mf.Vm[ir,1]=0.0
        mf.DeltaPot[ir,1]=12.0/sqrt(nucl.A)

        if r > nucl.R
            mf.Vc[ir]=(nucl.Z-1)/r*hb*fs
        else
            mf.Vc[ir]=(nucl.Z-1)*hb*fs*(3nucl.R^2-r^2)/2/nucl.R^3
        end
    end
end

function Vcent(l, r, nucl::Nucleus)
    return l*(l+1)*hb2/2/nucl.nuclmassc/r^2
end

function CalcPot(mf::MeanField, nucl::Nucleus, j, l, iq)
    V=zeros(Float64, ngrid)
    xls=(j/2)*(j/2+1)-l*(l+1)-0.75
    for ir=1:ngrid
        V[ir]=(mf.Vn[ir,iq]+mf.Vm[ir,iq]+mf.Vls[ir,iq]*xls
        +(mf.Vc[ir]+mf.Vcex[ir])*(2-iq)
        +Vcent(l,ir*dr, nucl)*mf.Binv[ir,iq]/nucl.Binv0)
    end
    return V
end

function Numerov(it, mf::MeanField, nucl::Nucleus, V, e, l, iq)
    phi=zeros(Float64, ngrid)
    fac=dr^2*(2*nucl.nuclmassc/hb2)

    phi0=dr^(l+1)
    phi1=(2+5.0/6.0*fac*nucl.Binv0/mf.Binv[1,iq]*(V[1]-e))*phi0
    phi1=phi1/(1-fac/12*nucl.Binv0/mf.Binv[2,iq]*(V[2]-e))

    phi[1]=phi0
    phi[2]=phi1

    for ir=3:ngrid
        dd=phi0-2*phi1
        dd-=fac/12*nucl.Binv0/mf.Binv[ir-2,iq]*(V[ir-2]-e)*phi0
        dd-=5*fac/6*nucl.Binv0/mf.Binv[ir-1,iq]*(V[ir-1]-e)*phi1

        cc=1.0-fac/12*nucl.Binv0/mf.Binv[ir,iq]*(V[ir]-e)
        cin=1.0/cc

        phi[ir]=-cin*dd

        phi0=phi1
        phi1=phi[ir]
    end

    #Normalization
    norm=0.0
    for ir=1:ngrid
        norm+=phi[ir]^2*dr*nucl.Binv0/mf.Binv[ir,iq]
    end

    nod=0
    for ir=1:ngrid
        phi[ir]=phi[ir]/sqrt(norm)*sqrt(nucl.Binv0/mf.Binv[ir,iq])
        if ir > 5 && phi[ir]*phi[ir-1] < 0
            nod+=1
        end
    end
    return nod, phi
end

function Numerov2(mf::MeanField, nucl::Nucleus, V, e, psi0, iq)
    phi=zeros(Float64, ngrid)
    fac=dr^2*(2*nucl.nuclmassc/hb2)
    ak=sqrt(2*nucl.nuclmassc/hb2*abs(e))

    phi[ngrid]   = exp(-ak*(ngrid  )*dr)
    phi[ngrid-1] = exp(-ak*(ngrid-1)*dr)

    Cr=1.5*nucl.A^(1.0/3.0)
    irc=round(Int64, Cr/dr)

    phi0=phi[ngrid]
    phi1=phi[ngrid-1]
    for ir=ngrid-2:-1:irc
        dd=phi0-2*phi1
        dd=dd-fac/12*nucl.Binv0/mf.Binv[ir+2,iq]*(V[ir+2]-e)*phi0
        dd=dd-5*fac/6*nucl.Binv0/mf.Binv[ir+1,iq]*(V[ir+1]-e)*phi1

        cc=-fac/12*nucl.Binv0/mf.Binv[ir,iq]*(V[ir]-e)+1.0
        cin=1.0/cc

        phi[ir]=-cin*dd

        phi0=phi1
        phi1=phi[ir]
    end

    #matching
    b =sqrt(nucl.Binv0/mf.Binv[irc,iq])
    bh=sqrt(nucl.Binv0/mf.Binv[irc+1,iq])
    ac=log((psi0[irc+1]*phi[irc]*b)/(psi0[irc]*phi[irc+1]*bh))/dr

    cmatch=sqrt(mf.Binv[irc,iq]/nucl.Binv0)*psi0[irc]/phi[irc]
    for ir=irc:ngrid
        psi0[ir]=phi[ir]*cmatch*sqrt(nucl.Binv0/mf.Binv[ir,iq])
    end

    #normalization
    norm=0.0
    for ir=1:ngrid
        norm+=psi0[ir]^2*dr
    end
    for ir=1:ngrid
        psi0[ir]/=sqrt(norm)
    end
    return psi0
end

function Spbasis(nucl::Nucleus, mf::MeanField)
    tol=10^(-6)
    V=zeros(Float64, ngrid)
    psi0=zeros(Float64, ngrid)
    for iq=2:-1:1
        nstate=0
        for l=0:lcut
            for j = abs(2*l-1) : 2: 2*l+1
                V=CalcPot(mf,nucl,j,l,iq)
                e=ecut
                nodemax0, psi0=Numerov(0, mf, nucl, V, e, l, iq)
                nodemax0-=1
                if nodemax0 > nodemax
                    nodemax0=nodemax
                end

                emin=-100.0
                emax=ecut
                for node0=0:nodemax0
                    nn=50
                    for i0=1:nn
                        e=(emin+emax)/2
                        node, psi0=Numerov(i0, mf, nucl, V, e, l, iq)
                        if node > node0
                            emax=e
                        else
                            emin=e
                        end
                    end

                    nstate+=1
                    if nstate > nstatemax
                        println("Increase nstatemax!")
                        exit()
                    end

                    if e < 0
                        psi0=Numerov2(mf, nucl, V, e, psi0, iq)
                    end

                    nucl.eh[nstate,iq]=e
                    nucl.lh[nstate,iq]=l
                    nucl.jh[nstate,iq]=j
                    nucl.node[nstate,iq]=node0
                    for ir=1:ngrid
                        nucl.psi[ir,nstate,iq]=psi0[ir]
                    end

                    emin=e
                    emax=ecut
                end
            end
        end
        nucl.nmax[iq]=nstate
    end
    sort_hf(nucl)
end

function Aux(nucl::Nucleus)
    df=zeros(Float64, ngrid, 2)
    ddf=zeros(Float64, ngrid, 2)
    dsj=zeros(Float64, ngrid, 2)
    drho=zeros(Float64, ngrid, 2)
    ddrho=zeros(Float64, ngrid, 2)
    for iq=1:2
        ddrho[:,iq].=Deriv25(nucl.rho[:,iq])
        for i=1:nucl.nmax[iq]
            psi=nucl.psi[:,i,iq]
            dpsi=Deriv5(psi)
            v2=nucl.v2[i,iq]
            j=nucl.jh[i,iq]
            l=nucl.lh[i,iq]
            for ir=1:ngrid # nabla^2 rho = rho''+2*rho'/2
                r=rgrid[ir]
                drho[ir,iq]+=v2*(j+1)*pi4m*(2*psi[ir]*dpsi[ir]/r^2-2*psi[ir]^2/r^3)
                dsj[ir,iq] +=v2*(j+1)*pi4m*(j/2*(j/2+1)-l*(l+1)-0.75)*(-3*psi[ir]^2/r^4+2*psi[ir]*dpsi[ir]/r^3)
            end
        end
        #
        for ir=1:ngrid
            r=rgrid[ir]
            df[ir,iq]=dsj[ir,iq]+2.0/r*nucl.sj[ir,iq] # nabla cdot sj
            ddf[ir,iq]=ddrho[ir,iq]+2.0/r*drho[ir,iq]
        end
    end
    return drho, ddrho, dsj, df, ddf
end

function SkyrmePot(nucl::Nucleus, skm::Skyrme)
    t0=skm.t0
    t1=skm.t1
    t2=skm.t2
    t3=skm.t3
    x0=skm.x0
    x1=skm.x1
    x2=skm.x2
    x3=skm.x3
    w0=skm.w0
    w0p=skm.w0p
    alp=skm.alp
#
    Vn=zeros(Float64, ngrid, 2)
    Vls=zeros(Float64, ngrid, 2)
    Vm=zeros(Float64, ngrid, 2)
    DeltaPot=zeros(Float64, ngrid, 2)
    f1=zeros(Float64, ngrid)
    f2=zeros(Float64, ngrid)
    f3=zeros(Float64, ngrid)
    Binv=zeros(Float64, ngrid, 2)
    dBinv=zeros(Float64, ngrid, 2)
    ddBinv=zeros(Float64, ngrid, 2)
    #
    drho, ddrho, dsj, df, ddf=Aux(nucl)
    #
    for iq=1:2
        for ir=1:ngrid # nabla^2 rho = rho''+2*rho'/2
            r=rgrid[ir]
            rhop=nucl.rho[ir,1]
            rhon=nucl.rho[ir,2]
            rho=nucl.rho[ir,1]+nucl.rho[ir,2]
            drhot=drho[ir,1]+drho[ir,2]
            tau=nucl.tau[ir,1]+nucl.tau[ir,2]
                                
            Vn[ir,iq]=(t0*(1+x0/2)*rho-t0*(0.5+x0)*nucl.rho[ir,iq]
                    +t3/12*(1+x3/2)*(alp+2)*rho^(alp+1)
                    -t3/12*(0.5+x3)*(alp*rho^(alp-1)*(rhon^2+rhop^2)+2*rho^alp*nucl.rho[ir,iq])
                    +(t1*(1+x1/2)+t2*(1+x2/2))*tau/4-(t1*(0.5+x1)-t2*(0.5+x2))*nucl.tau[ir,iq]/4
                    -(3*t1*(1+x1/2)-t2*(1+x2/2))*(ddf[ir,1]+ddf[ir,2])/8
                    +(3*t1*(0.5+x1)+t2*(0.5+x2))*ddf[ir,iq]/8
                    -0.5*w0*(df[ir,1]+df[ir,2])-0.5*w0p*df[ir,iq])
            Vls[ir,iq]=0.5*w0*drhot+0.5*w0p*drho[ir,iq]
            if skm.ij2
                Vls[ir,iq]+=((t1-t2)-(t1*x1+t2*x2))/8*sj[ir,iq]+(-(t1*x1+t2*x2))/8*sj[ir,3-iq]
            end
            Vls[ir,iq]=Vls[ir,iq]/r

            Binv[ir,iq]=(hb2/2/nucl.nuclmassc
               +0.25*(t1*(1+x1/2)+t2*(1+x2/2))*rho
               -0.25*(t1*(0.5+x1)-t2*(0.5+x2))*nucl.rho[ir,iq])
          
            dBinv[ir,iq] =(
          (t1*(1+x1/2)+t2*(1+x2/2))*drhot/4
         -(t1*(0.5+x1)-t2*(0.5+x2))*drho[ir,iq]/4)
            ddBinv[ir,iq]=(t1*(1+x1/2)+t2*(1+x2/2))*(ddrho[ir,1]+ddrho[ir,2])/4-(t1*(0.5+x1)-t2*(0.5+x2))*ddrho[ir,iq]/4

           DeltaPot[ir,iq]=skm.vpair[iq]/2*nucl.rhopair[ir,iq]
        end
    end
#     potential due to the effective mass
    for iq=1:2
        for ir=1:ngrid
            r=rgrid[ir]

            f1[ir]=Binv[ir,iq]/nucl.Binv0
            f2[ir]=-sqrt(nucl.Binv0)/2*Binv[ir,iq]^(-1.5)*dBinv[ir,iq] 
            # d/dr sqrt(binv0/binv(iq,ir))
            f3[ir]=3*sqrt(nucl.Binv0)/4*Binv[ir,iq]^(-2.5)*dBinv[ir,iq]^2-sqrt(nucl.Binv0)/2*Binv[ir,iq]^(-1.5)*ddBinv[ir,iq] 
            # d^2/dr^2 sqrt(binv0/binv(iq,ir))
            Vm[ir,iq]=-nucl.Binv0*sqrt(Binv[ir,iq]/nucl.Binv0)*(dBinv[ir,iq]/nucl.Binv0*f2[ir]+f1[ir]*(f3[ir]+2.0/r*f2[ir]))
        end
    end
    return Vn, Vls, Vm, Binv, DeltaPot
end

function Potential(mf::MeanField, nucl::Nucleus, skm::Skyrme, xmu, it)
    if it == 1
        WoodSaxon(mf, nucl)
    else
        Vn,Vls,Vm,Binv,DeltaPot=SkyrmePot(nucl, skm)
        Vc,Vcex=Coulomb(nucl, mf)
        for iq=1:2
            for ir=1:ngrid
                mf.Vn[ir,iq]  = (1-xmu)*Vn[ir,iq]  + xmu*mf.Vn[ir,iq]
                mf.Vls[ir,iq] = (1-xmu)*Vls[ir,iq] + xmu*mf.Vls[ir,iq]
                mf.Vm[ir,iq]  = Vm[ir,iq]
                mf.Binv[ir,iq]= Binv[ir,iq]
                mf.DeltaPot[ir,iq]= DeltaPot[ir,iq]
                mf.Vc[ir]  = Vc[ir]
                mf.Vcex[ir]= Vcex[ir]
            end
        end
    end
end

function sort_hf(s::Nucleus)
    for iq=1:2
        nmax=s.nmax[iq]

        for i=1:nmax-1
            kk=i
            p=s.eh[i,iq]
            for j=i+1:nmax
                if s.eh[j,iq] <= p
                    kk=j
                    p=s.eh[j,iq]
                end
            end

            if kk != i
                s.eh[kk,iq]=s.eh[i,iq]
                s.eh[i,iq]=p

                p=s.jh[i,iq]
                s.jh[i,iq]=s.jh[kk,iq]
                s.jh[kk,iq]=p

                p=s.v2[i,iq]
                s.v2[i,iq]=s.v2[kk,iq]
                s.v2[kk,iq]=p

                p=s.lh[i,iq]
                s.lh[i,iq]=s.lh[kk,iq]
                s.lh[kk,iq]=p

                p=s.node[i,iq]
                s.node[i,iq]=s.node[kk,iq]
                s.node[kk,iq]=p

                for ir=1:ngrid
                    p=s.psi[ir,i,iq]
                    s.psi[ir,i,iq]=s.psi[ir,kk,iq]
                    s.psi[ir,kk,iq]=p
                end
            end
        end
    end
end

function gap(nucl::Nucleus, mf::MeanField, iq)
    nucl.delta=zeros(Float64,ngrid,2)
    for i=1:nucl.nmax[iq]
        for ir=1:ngrid
            nucl.delta[i,iq]+=mf.DeltaPot[ir,iq]*nucl.psi[ir,i,iq]^2*dr
        end
    end
end

function bcs(nucl::Nucleus, mf::MeanField)
    tol=1e-6
    for iq=1:2
        if iq == 1
            N=nucl.Z
        else
            N=nucl.N
        end

        gap(nucl, mf, iq)

        efmax = 20.0
        efmin =-60.0
	eftry =  0.0
        for ite=1:100
            eftry=(efmax+efmin)/2

            sum=0.0
            for i=1:nucl.nmax[iq]
                eqp=sqrt((nucl.eh[i,iq]-eftry)^2+(nucl.delta[i,iq])^2)
                nucl.v2[i,iq]=0.5*(1-(nucl.eh[i,iq]-eftry)/eqp)
                sum+=(nucl.jh[i,iq]+1)*nucl.v2[i,iq]
            end

            if abs(sum-N)/N < tol
                break
            elseif sum > N
                efmax=eftry
            elseif sum < N
                efmin=eftry
            end
        end
        nucl.eferm[iq]=eftry
    end
end

function Coulomb(nucl::Nucleus, mf::MeanField)
    Vc=zeros(Float64,ngrid)
    Vcex=zeros(Float64,ngrid)
    s=zeros(Float64,ngrid)
    for ir=1:ngrid
        s[ir]=-4pi*hb*fs*nucl.rho[ir,1]*rgrid[ir]
    end
    #
    charge=0.0
    for ir=1:ngrid
        r=rgrid[ir]
        charge+=4pi*r^2*hb*fs*nucl.rho[ir,1]*dr
    end
    #
    Vc0=zeros(Float64,ngrid)
    Vc0[ngrid]=charge
    Vc0[ngrid-1]=charge
    for ir=ngrid-2:-1:1
        Vc0[ir]=2*Vc0[ir+1]-Vc0[ir+2]+dr^2/12*(s[ir]+10*s[ir+1]+s[ir+2])
    end
    #
    for ir=1:ngrid
        Vc[ir]=Vc0[ir]/rgrid[ir]
    end
    #exchange term
    for ir=1:ngrid
        Vcex[ir]=-(3.0/pi)^(1.0/3.0)*hb*fs*nucl.rho[ir,1]^(1.0/3.0)
    end
    return Vc, Vcex
end

function density(nucl::Nucleus)
    for iq=1:2
        nucl.rho[:,iq].=0.0
        nucl.tau[:,iq].=0.0
        nucl.sj[:,iq].=0.0
        nucl.rhopair[:,iq].=0.0
        for i=1:nucl.nmax[iq]
            j=nucl.jh[i,iq]
            l=nucl.lh[i,iq]
            v2=nucl.v2[i,iq]
            dpsi=Deriv5(nucl.psi[:,i,iq])
            nucl.eqp[i,iq]=sqrt((nucl.eh[i,iq]-nucl.eferm[iq])^2+(nucl.delta[i,iq])^2)
            for ir=1:ngrid
                r=rgrid[ir]
                r2=r^2
                psi=nucl.psi[ir,i,iq]
                nucl.rho[ir,iq]+=v2*(j+1)*pi4m*(psi/r)^2
                nucl.tau[ir,iq]+=v2*(j+1)*pi4m*((dpsi[ir]/r-psi/r2)^2+l*(l+1)/r2*(psi/r)^2)
                nucl.sj[ir,iq] +=v2*(j+1)*pi4m*(j/2*(j/2+1)-l*(l+1)-0.75)/r*(psi/r)^2
                nucl.rhopair[ir,iq]+=-sqrt(v2*(1.0-v2))*(j+1)*pi4m*(psi/r)^2
            end
        end
    end
    for iq=1:2
        if iq == 1
            N=nucl.Z
        else
            N=nucl.N
        end
        rms=0.0
        for ir=1:ngrid
            rms+=4π*nucl.rho[ir,iq]*(ir*dr)^4*dr
        end
        nucl.radii[iq]=sqrt(rms/N)
    end
    nucl.radii[3]=sqrt((nucl.radii[1]^2*nucl.Z+nucl.radii[2]^2*nucl.N)/(nucl.A))
end

function energy(mf::MeanField, nucl::Nucleus, skm::Skyrme)
    #Single Particle Energy
    esp=0.0
    for iq=1:2
        for i=1:nucl.nmax[iq]
            esp+=nucl.v2[i,iq]*(nucl.jh[i,iq]+1)*nucl.eh[i,iq]
        end
    end
    #Kinetic Energy
    ekin=0.0
    for ir=1:ngrid
        r=rgrid[ir]
        ekin+=4pi*r^2*hb2/2/nucl.nuclmassc*(nucl.tau[ir,1]+nucl.tau[ir,2])*dr
    end
    #Rearrangement term
    erearrangement=0.0
    for ir=1:ngrid
        r=rgrid[ir]
        r2=r^2
        rho=nucl.rho[ir,1]+nucl.rho[ir,2]
        erearrangement-=4pi*r2*(
        skm.alp/24*skm.t3*rho^skm.alp*((1+skm.x3/2)*rho^2
        -(0.5+skm.x3)*(nucl.rho[ir,1]^2+nucl.rho[ir,2]^2)))*dr
        erearrangement-=4pi*r2*(0.25*(3/pi)^(1.0/3.0)*hb*fs*nucl.rho[ir,1]^(4.0/3.0))*dr
    end
    #Pairing Energy
    epair=zeros(Float64,2)
    for iq=1:2
        for ir=1:ngrid
            r=rgrid[ir]
            epair[iq]+=4pi*r^2*mf.DeltaPot[ir,iq]*nucl.rhopair[ir,iq]/2*dr
        end
    end
    nucl.Etot=(esp+ekin)/2+erearrangement+epair[1]+epair[2]
    nucl.Epair=epair
end

function Result(nucl::Nucleus, mf::MeanField)
    Write_levels(nucl)
    Write_potential(mf)
    Write_density(nucl)
end

function Write_potential(mf::MeanField)
    f=open("potential.out","w")
    println(f,"      r   proton Vn   neutron Vn")
    for ir=1:ngrid
        r=ir*dr
        @printf(f,"%6.3f %12.7f %12.7f\n",r,mf.Vn[ir,1],mf.Vn[ir,2])
    end
    close(f)
end

function Derivatives(nucl::Nucleus)
    e=zeros(Float64,ngrid)
    for iq=1:2
        for i=1:nucl.nmax[iq]
            for ir=1:ngrid
                e[ir]=nucl.psi[ir,i,iq]
            end
            nucl.dpsi[:,i,iq].=Deriv5(e)
            nucl.ddpsi[:,i,iq].=Deriv25(e)
        end
    end
    return
end

function Calc_shfbcs(nucl::Nucleus, skm::Skyrme, mf::MeanField)
    itmax=100
    xmu=0.9
    if iprint
        @printf("Iteration starts\n")
        @printf("  Times    Etot  Entropy    efermn    efermp    epairn    epairp\n")
    end
    for it=1:itmax
        if it == 15
            xmu=0.9
        elseif it == 20
            xmu=0.9
        elseif it == 25
            xmu=0.9
        elseif it == 30
            xmu=0.9
        end
        E0=nucl.Etot
        Potential(mf, nucl, skm, xmu, it)
        Spbasis(nucl, mf)
        bcs(nucl, mf)
        density(nucl)
        energy(mf, nucl, skm)
        acc=(nucl.Etot-E0)/E0
    if iprint
            @printf("%5.0f%10.4f%9.4f%10.4f%10.4f%10.4f%10.4f%12.4e%9.3f\n",
                     it,nucl.Etot,nucl.Entropy,nucl.eferm[2],nucl.eferm[1],
                        nucl.Epair[2],nucl.Epair[1],acc,xmu)
    end
        if it > 15 && abs(acc) < tolk
            break
        end
    end
    Result(nucl,mf)
    Derivatives(nucl)
end

end #module shfbcs

