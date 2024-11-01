module Angle

using CGcoefficient
using ..Basic

export Yl, Yjl, YLSjl, Ynlsig

function Yl(L,l1,l2)
    if l1 < 0 || l2 < 0
        return 0
    else
        return ((-1)^(l1)*sqrt((2l1+1)*(2l2+1)*(2L+1)*pi4m)
                *float(threeJ(l1,L,l2,0,0,0)))
    end
end

function Yjl(L,j1,l1,j2,l2)
    if l1 < 0 || l2 < 0 || j1 > 2l1+1 || j1 < 2l1-1 || j2 > 2l2+1 || j2 < 2l2-1
        return 0
    end
    fac=(j1+1)*(j2+1)*(2l1+1)*(2l2+1)*(2L+1)*pi4m
    return ((-1)^((1+j2)/2+L)*sqrt(fac)
            *float(threeJ(l1, L, l2, 0, 0, 0))
            *float(sixJ(l1, j1//2, 1//2, j2//2, l2, L)))
end
#<J1 L1 || [Y_L S]^J || J2 L2>
function YLSjl(j1::Int, l1::Int, j2::Int, l2::Int, J::Int, L::Int)
    if L < 0 || l1 < 0 || l2 < 0 || j1 > 2l1+1 || j1 < 2l1-1 || j2 > 2l2+1 || j2 < 2l2-1
        return 0
    end
    fac=6*(j1+1)*(j2+1)*(2J+1)*(2l1+1)*(2l2+1)*(2L+1)*pi4m
    return (sqrt(fac)
            *float(nineJ(l1,l2,L,1//2,1//2,1,j1//2,j2//2,J))
            *(-1)^(l1)
            *float(threeJ(l1,L,l2,0,0,0)))
end

function Ynlsig(JJ::Int, J1::Int, L1::Int, J2::Int, L2::Int, sign::Int)
    fac=(-1)^(L1+L2)*sqrt(6.0*(J1+1)*(J2+1)*(2*JJ+1))
    #LL=JJ+1
    Ynlsig=fac*(-1)^(JJ+1)*(sqrt(2*(JJ+1)+1)
                 *float(nineJ(L1,L2,JJ+1,1//2,1//2,1,J1//2,J2//2,JJ))
                 *float(sixJ(JJ,1,JJ+1,L2,L1,L2+sign*1))*Yl(JJ,L1,L2+sign*1))
    #LL=JJ-1
    if JJ-1 >= 0
        Ynlsig+=fac*(-1)^(JJ-1)*(sqrt(2*(JJ-1)+1)
                     *float(nineJ(L1,L2,JJ-1,1//2,1//2,1,J1//2,J2//2,JJ))
                     *float(sixJ(JJ,1,JJ-1,L2,L1,L2+sign*1))*Yl(JJ,L1,L2+sign*1))
    end
    return Ynlsig
end

# < J1 L1 || [[Y_LL nab]^JJ || J2 L2 >
function Ynab(LL,J1,L1,J2,L2,JJ)
    Ynab=zeros(Float64,2)
    if LL < 0
        return ynab
    end
    fac=(-1)^((1+J2)/2+L2)*(sqrt((J1+1.0)*(J2+1.0)*(2*JJ+1.0)) 
                            *float(sixj(L1,J1//2,1//2,J2//2,L2,JJ)))
    Ynab[1]=fac*r6j(LL,1,JJ,L2,L1,(L2-1))*Yl(LL,L1,L2-1)
    Ynab[2]=fac*r6j(LL,1,JJ,L2,L1,(L2+1))*Yl(LL,L1,L2+1)
    return ynab
end

end
