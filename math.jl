module Math
export Deriv5, Deriv25, Deriv5i
using ..Basic

function Deriv5(f)
    df=zeros(Float64, ngrid)
    dr120=1.0/(120*dr)
    df[1]=(-274*f[1]+600*f[2]-600*f[3]+400*f[4]-150*f[5]+ 24*f[6])*dr120
    df[2]=(- 24*f[1]-130*f[2]+240*f[3]-120*f[4]+ 40*f[5]-  6*f[6])*dr120
    df[3]=(   6*f[1]- 60*f[2]- 40*f[3]+120*f[4]- 30*f[5]+  4*f[6])*dr120
    df[4]=(   6*f[2]- 60*f[3]- 40*f[4]+120*f[5]- 30*f[6]+  4*f[7])*dr120

    df[ngrid-3]=( -4*f[ngrid-6]+ 30*f[ngrid-5]-120*f[ngrid-4]
                 +40*f[ngrid-3]+ 60*f[ngrid-2]-  6*f[ngrid-1])*dr120
    df[ngrid-2]=( -4*f[ngrid-5]+ 30*f[ngrid-4]-120*f[ngrid-3]
                 +40*f[ngrid-2]+ 60*f[ngrid-1]-  6*f[ngrid  ])*dr120
    df[ngrid-1]=(-24*f[ngrid-5]+150*f[ngrid-4]-400*f[ngrid-3]
                +600*f[ngrid-2]-600*f[ngrid-1]+274*f[ngrid  ])*dr120
    df[ngrid  ]=(  6*f[ngrid-5]- 40*f[ngrid-4]+120*f[ngrid-3]
                -240*f[ngrid-2]+130*f[ngrid-1]+ 24*f[ngrid  ])*dr120
    dr840=1.0/(840*dr)
    for ir=5:ngrid-4
        df[ir]=( -3*(f[ir+4]-f[ir-4]) +32*(f[ir+3]-f[ir-3])
               -168*(f[ir+2]-f[ir-2])+672*(f[ir+1]-f[ir-1]))*dr840
    end
    return df
end

function Deriv25(f)
    ddf=zeros(Float64, ngrid)
    dr2=dr^2
    dr260=1.0/(60*dr2)
    
    ddf[1]=(225*f[1]-770*f[2]+1070*f[3]-780*f[4]+305*f[5]- 50*f[6])*dr260
    ddf[2]=( 50*f[1]- 75*f[2]-  20*f[3]+ 70*f[4]- 30*f[5]+  5*f[6])*dr260
    ddf[3]=(- 5*f[1]+ 80*f[2]- 150*f[3]+ 80*f[4]-  5*f[5])*dr260
    ddf[4]=(- 5*f[2]+ 80*f[3]- 150*f[4]+ 80*f[5]-  5*f[6])*dr260

    ddf[ngrid-3]=(- 5*f[ngrid-5]+ 80*f[ngrid-4]-150*f[ngrid-3]+  80*f[ngrid-2]-  5*f[ngrid-1])*dr260
    ddf[ngrid-2]=(- 5*f[ngrid-4]+ 80*f[ngrid-3]-150*f[ngrid-2]+  80*f[ngrid-1]-  5*f[ngrid  ])*dr260
    ddf[ngrid-1]=(  5*f[ngrid-5]- 30*f[ngrid-4]+ 70*f[ngrid-3]-  20*f[ngrid-2]- 75*f[ngrid-1]+ 50*f[ngrid])*dr260
    ddf[ngrid  ]=(-50*f[ngrid-5]+305*f[ngrid-4]-780*f[ngrid-3]+1070*f[ngrid-2]-770*f[ngrid-1]+225*f[ngrid])*dr260

    dr212=1.0/(12*dr2)
    for ir=3:ngrid-4
        ddf[ir]=-f[ir-2]+16*f[ir-1]-30*f[ir]+16*f[ir+1]-f[ir+2]
        ddf[ir]=ddf[ir]*dr212
    end
    return ddf
end

function Deriv5i(f)
    df=zeros(Complex{Float64}, ngrid)
    dr120=1.0/(120*dr)
    df[1]=(-274*f[1]+600*f[2]-600*f[3]+400*f[4]-150*f[5]+ 24*f[6])*dr120
    df[2]=(- 24*f[1]-130*f[2]+240*f[3]-120*f[4]+ 40*f[5]-  6*f[6])*dr120
    df[3]=(   6*f[1]- 60*f[2]- 40*f[3]+120*f[4]- 30*f[5]+  4*f[6])*dr120
    df[4]=(   6*f[2]- 60*f[3]- 40*f[4]+120*f[5]- 30*f[6]+  4*f[7])*dr120

    df[ngrid-3]=( -4*f[ngrid-6]+ 30*f[ngrid-5]-120*f[ngrid-4]
                 +40*f[ngrid-3]+ 60*f[ngrid-2]-  6*f[ngrid-1])*dr120
    df[ngrid-2]=( -4*f[ngrid-5]+ 30*f[ngrid-4]-120*f[ngrid-3]
                 +40*f[ngrid-2]+ 60*f[ngrid-1]-  6*f[ngrid  ])*dr120
    df[ngrid-1]=(-24*f[ngrid-5]+150*f[ngrid-4]-400*f[ngrid-3]
                +600*f[ngrid-2]-600*f[ngrid-1]+274*f[ngrid  ])*dr120
    df[ngrid  ]=(  6*f[ngrid-5]- 40*f[ngrid-4]+120*f[ngrid-3]
                -240*f[ngrid-2]+130*f[ngrid-1]+ 24*f[ngrid  ])*dr120
    dr840=1.0/(840*dr)
    for ir=5:ngrid-4
        df[ir]=( -3*(f[ir+4]-f[ir-4]) +32*(f[ir+3]-f[ir-3])
               -168*(f[ir+2]-f[ir-2])+672*(f[ir+1]-f[ir-1]))*dr840
    end
    return df
end

end