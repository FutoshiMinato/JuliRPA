module Basic
export hb,hb2,hb22m,pi4m,fs,nuclmass,nstatemax,lcut,nodemax,ecut,rbox,dr,ngrid,tolk
export fs
const hb=197.3285879
const hb2=hb*hb
const pi4m=1.0/(4*pi)
const fs=1.0/137.03604
const hb22m=20.73552985
const nuclmass=hb2/2/hb22m

const nstatemax=20
const lcut=5
const nodemax=3
const ecut=0.0
const rbox=12
const dr=0.1
const ngrid=Int(rbox/dr)
const tolk=1e-5

const fs=1.0/(137.03604)

end
