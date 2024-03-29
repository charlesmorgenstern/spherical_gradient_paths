using Plots; plotlyjs() #good for 3d plots
#using Plots             #good for 2d plots
using LinearAlgebra
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
function plotfunc3d()

f(t,p)=cos.(p).*sin.(t)+cos.(2 .*t )

n = 300
t = LinRange(0,pi,n)
p = LinRange(0,2*pi,n)

tt = t*ones(n)'
pp = ones(n)*p'

xx = @. sin(t)*cos(p')
yy = @. sin(t)*sin(p')
zz = @. cos(t)*ones(n)'

fig=surface(xx, yy, zz, fill_z = f(tt,pp) , size=(600,600))

return fig
end
#############################################################
#############################################################
function plotfunc2d()

f(t,p)=cos(p)*sin(t)+cos(2*t)

p = range(0, 2*pi, length=100)
t = range(0, pi, length=100)

z = @. f(t', p)

plt=contour(t, p, z, levels=40, color=:plasma, xlabel="theta", ylabel="phi")
#plt=heatmap(t, p, z, color=:plasma, xlabel="phi", ylabel="theta")

return plt
end
#############################################################
#############################################################
function plotfunc2ddom()

f(t,p)=cos(p)*sin(t)+cos(2*t)

p = range(0-pi/2, 2*pi+pi/2, length=100)
t = range(0-pi/2, pi+pi/2, length=100)

z = @. f(t', p)

plt=contour(t, p, z, levels=40, color=:plasma, xlabel="theta", ylabel="phi")
#plt=heatmap(t, p, z, color=:plasma, xlabel="phi", ylabel="theta")

plot!([0, pi, pi, 0, 0],[0, 0, 2*pi, 2*pi, 0],legend=false)

return plt
end
#############################################################
#############################################################
function plotcritpts2d()

f(t,p)=cos(p)*sin(t)+cos(2*t)

p = range(0, 2*pi, length=100)
t = range(0, pi, length=100)

z = @. f(t', p)

pts=[pi/2 pi;
     pi/2 0;
     pi/2 2*pi;
     2*atan(4+sqrt(15))   0;
     2*atan(4-sqrt(15))   0;
     2*atan(4+sqrt(15))   2*pi;
     2*atan(4-sqrt(15))   2*pi]


plt=contour(t, p, z, levels=40, color=:plasma, xlabel="theta", ylabel="phi")
#plt=heatmap(t, p, z, color=:plasma, xlabel="phi", ylabel="theta")
scatter!(pts[:,1],pts[:,2],legend=false)

return plt
end
#############################################################
#############################################################
function grad(t,p)

g=[cos(p)*cos(t)-2*sin(2*t),-sin(p)]

return g
end
#############################################################
#############################################################
function gradpath(t,p)

h=.01

cps=[pi/2 pi;
     pi/2 0;
     pi/2 2*pi;
     2*atan(4+sqrt(15))   0;
     2*atan(4-sqrt(15))   0;
     2*atan(4+sqrt(15))   2*pi;
     2*atan(4-sqrt(15))   2*pi]


gp1=[t p]
gp2=[t p]
check=true
tn=t
pn=p
cutoff=1000

while check
g=-grad(tn,pn)
g=g./norm(g)
tt=tn+h*g[1]
pp=pn+h*g[2]

if tt<0
tt=0
pp=pp+pi
end
if tt>pi
tt=pi
pp=pp+pi
end

ncp=length(cps[:,1])
for j=1:ncp
if norm([tt pp]-cps[j,:]')<h
check=false
tt=cps[j,1]
pp=cps[j,2]
break
end
end

gp1=[gp1;tt pp]
tn=tt
pn=pp
end

tn=t
pn=p
check=true
while check
g=-grad(tn,pn)
g=g./norm(g)
tt=tn-h*g[1]
pp=pn-h*g[2]

if tt<0
tt=0
pp=pp+pi
end
if tt>pi
tt=pi
pp=pp+pi
end

ncp=length(cps[:,1])
for j=1:ncp
if norm([tt pp]-cps[j,:]')<h
check=false
tt=cps[j,1]
pp=cps[j,2]
break
end
end

gp2=[gp2;tt pp]
tn=tt
pn=pp
end

gp1=gp1[2:end,:]
gp1=reverse(gp1,dims=1)

gp=[gp1;gp2]

return gp
end
#############################################################
#############################################################
function gradpathrk4(t,p)

h=.01

cps=[pi/2 pi;
     pi/2 0;
     pi/2 2*pi;
     2*atan(4+sqrt(15))   0;
     2*atan(4-sqrt(15))   0;
     2*atan(4+sqrt(15))   2*pi;
     2*atan(4-sqrt(15))   2*pi]


gp1=[t p]
gp2=[t p]
check=true
tn=t
pn=p
cutoff=1000

while check
g=-grad(tn,pn)
g=g./norm(g)

 kt1=g[1]
 kp1=g[2]
 g=-grad(tn+h*kt1/2,pn+h*kp1/2)
 g=g./norm(g)
 kt2=g[1]
 kp2=g[2]
 g=-grad(tn+h*kt2/2,pn+h*kp2/2)
 g=g./norm(g)
 kt3=g[1]
 kp3=g[2]
 g=-grad(tn+h*kt3,pn+h*kp3)
 g=g./norm(g)
 kt4=g[1]
 kp4=g[2]


 tt=tn+h/6*(kt1+2*kt2+2*kt3+kt4)
 pp=pn+h/6*(kp1+2*kp2+2*kp3+kp4)

if tt<0
tt=0
pp=pp+pi
end
if tt>pi
tt=pi
pp=pp+pi
end

ncp=length(cps[:,1])
for j=1:ncp
if norm([tt pp]-cps[j,:]')<h
check=false
tt=cps[j,1]
pp=cps[j,2]
break
end
end

gp1=[gp1;tt pp]
tn=tt
pn=pp
end

tn=t
pn=p
check=true
while check
g=grad(tn,pn)
g=g./norm(g)

 kt1=g[1]
 kp1=g[2]
 g=grad(tn+h*kt1/2,pn+h*kp1/2)
 g=g./norm(g)
 kt2=g[1]
 kp2=g[2]
 g=grad(tn+h*kt2/2,pn+h*kp2/2)
 g=g./norm(g)
 kt3=g[1]
 kp3=g[2]
 g=grad(tn+h*kt3,pn+h*kp3)
 g=g./norm(g)
 kt4=g[1]
 kp4=g[2]


 tt=tn+h/6*(kt1+2*kt2+2*kt3+kt4)
 pp=pn+h/6*(kp1+2*kp2+2*kp3+kp4)


if tt<0
tt=0
pp=pp+pi
end
if tt>pi
tt=pi
pp=pp+pi
end

ncp=length(cps[:,1])
for j=1:ncp
if norm([tt pp]-cps[j,:]')<h
check=false
tt=cps[j,1]
pp=cps[j,2]
break
end
end

gp2=[gp2;tt pp]
tn=tt
pn=pp
end

gp1=gp1[2:end,:]
gp1=reverse(gp1,dims=1)

gp=[gp1;gp2]

return gp
end
#############################################################
#############################################################
function spheretocart(gp)
np=length(gp[:,1])
out=zeros(np,3)

r=1.01
for i=1:np
t=gp[i,1]
p=gp[i,2]
out[i,:]=[r*sin(t)*cos(p) r*sin(t)*sin(p) r*cos(t)]
end

return out
end
#############################################################
#############################################################
function getphis(t)

f(t,p)=cos.(p).*sin.(t)+cos.(2 .*t)

n=10000
h=.0001

p=collect(range(0,2*pi,n))
p=p[1:end-1]
v=zeros(n-1)

for i=1:n-1
if t==0
v[i]=f(t+h,p[i])
else
v[i]=f(t-h,p[i])
end
end

d1=findmin(v)[2]
d2=findmax(v)[2]
p1=p[d1]
p2=p[d2]

out=[p1 p2]
return out
end
#############################################################
#############################################################
function plotpath(t1,p1)
f(t,p)=cos(p)*sin(t)+cos(2*t)

p = range(0, 2*pi, length=100)
t = range(0, pi, length=100)

z = @. f(t', p)

gp=gradpath(t1,p1)
fig=contour(t, p, z, levels=40, color=:plasma, xlabel="theta", ylabel="phi")
plot!(gp[:,1],gp[:,2],legend=false)
scatter!([t1],[p1],legend=false)

return fig
end
#############################################################
#############################################################
function plotpathrk4(t1,p1)
f(t,p)=cos(p)*sin(t)+cos(2*t)

p = range(0, 2*pi, length=100)
t = range(0, pi, length=100)

z = @. f(t', p)

gp=gradpathrk4(t1,p1)
fig=contour(t, p, z, levels=40, color=:plasma, xlabel="theta", ylabel="phi")
plot!(gp[:,1],gp[:,2],legend=false)
scatter!([t1],[p1],legend=false)

return fig
end
#############################################################
#############################################################
function plotpath3d(t1,p1)

f(t,p)=cos.(p).*sin.(t)+cos.(2 .*t)

n = 300
t = LinRange(0,pi,n)
p = LinRange(0,2*pi,n)

tt = t*ones(n)'
pp = ones(n)*p'

xx = @. sin(t)*cos(p')
yy = @. sin(t)*sin(p')
zz = @. cos(t)*ones(n)'

fig=surface(xx, yy, zz, fill_z = f(tt,pp) , size=(600,600))

gp=gradpath(t1,p1)
gp=spheretocart(gp)

plot!(gp[:,1], gp[:,2], gp[:,3], legend=false, linewidth=5)
scatter!([0 0],[0 0],[-1.01 1.01],legend=false)
return fig
end
#############################################################
#############################################################
function plotpath3drk4(t1,p1)

f(t,p)=cos.(p).*sin.(t)+cos.(2 .*t)

n = 300
t = LinRange(0,pi,n)
p = LinRange(0,2*pi,n)

tt = t*ones(n)'
pp = ones(n)*p'

xx = @. sin(t)*cos(p')
yy = @. sin(t)*sin(p')
zz = @. cos(t)*ones(n)'

fig=surface(xx, yy, zz, fill_z = f(tt,pp) , size=(600,600))

gp=gradpathrk4(t1,p1)
gp=spheretocart(gp)

plot!(gp[:,1], gp[:,2], gp[:,3], legend=false, linewidth=5)
scatter!([0 0],[0 0],[-1.01 1.01],legend=false)
return fig
end
#############################################################
#############################################################
function gradpathflag(t,p)
h=.001
flag=1

cps=[pi/2 pi;
     pi/2 0;
     pi/2 2*pi;
     2*atan(4+sqrt(15))   0;
     2*atan(4-sqrt(15))   0;
     2*atan(4+sqrt(15))   2*pi;
     2*atan(4-sqrt(15))   2*pi]

gp1=[t p]
gp2=[t p]
check=true
tn=t
pn=p

while check
if abs(tn)<10^-2 || abs(tn-pi)<10^-2
flag=0
end
g=-grad(tn,pn)
g=g./norm(g)
tt=tn+h*g[1]
pp=pn+h*g[2]

if tt<0
tt=0
pp=pp+pi
end
if tt>pi
tt=pi
pp=pp+pi
end

ncp=length(cps[:,1])
for j=1:ncp
if norm([tt pp]-cps[j,:]')<h
check=false
tt=cps[j,1]
pp=cps[j,2]
break
end
end

gp1=[gp1;tt pp]
tn=tt
pn=pp
end

tn=t
pn=p
check=true
cnt=1
while check
if abs(tn)<10^-2 || abs(tn-pi)<10^-2
flag=0
end
g=-grad(tn,pn)
g=g./norm(g)
tt=tn-h*g[1]
pp=pn-h*g[2]

if tt<0
tt=0
pp=pp+pi
end
if tt>pi
tt=pi
pp=pp+pi
end

ncp=length(cps[:,1])
for j=1:ncp
if norm([tt pp]-cps[j,:]')<h
check=false
tt=cps[j,1]
pp=cps[j,2]
break
end
end

gp2=[gp2;tt pp]
tn=tt
pn=pp
end

gp1=gp1[2:end,:]
gp1=reverse(gp1,dims=1)

gp=[gp1;gp2]

if flag==1
out=gp
else
out=0
end

return out
end
#############################################################
#############################################################
function gradpathflagrk4(t,p)
h=.001
flag=1

cps=[pi/2 pi;
     pi/2 0;
     pi/2 2*pi;
     2*atan(4+sqrt(15))   0;
     2*atan(4-sqrt(15))   0;
     2*atan(4+sqrt(15))   2*pi;
     2*atan(4-sqrt(15))   2*pi]

gp1=[t p]
gp2=[t p]
check=true
tn=t
pn=p

while check
if abs(tn)<10^-2 || abs(tn-pi)<10^-2
flag=0
end
g=-grad(tn,pn)
g=g./norm(g)

 kt1=g[1]
 kp1=g[2]
 g=-grad(tn+h*kt1/2,pn+h*kp1/2)
 g=g./norm(g)
 kt2=g[1]
 kp2=g[2]
 g=-grad(tn+h*kt2/2,pn+h*kp2/2)
 g=g./norm(g)
 kt3=g[1]
 kp3=g[2]
 g=-grad(tn+h*kt3,pn+h*kp3)
 g=g./norm(g)
 kt4=g[1]
 kp4=g[2]


 tt=tn+h/6*(kt1+2*kt2+2*kt3+kt4)
 pp=pn+h/6*(kp1+2*kp2+2*kp3+kp4)

if tt<0
tt=0
pp=pp+pi
end
if tt>pi
tt=pi
pp=pp+pi
end

ncp=length(cps[:,1])
for j=1:ncp
if norm([tt pp]-cps[j,:]')<h
check=false
tt=cps[j,1]
pp=cps[j,2]
break
end
end

gp1=[gp1;tt pp]
tn=tt
pn=pp
end

tn=t
pn=p
check=true
cnt=1
while check
if abs(tn)<10^-2 || abs(tn-pi)<10^-2
flag=0
end
g=grad(tn,pn)
g=g./norm(g)

 kt1=g[1]
 kp1=g[2]
 g=grad(tn+h*kt1/2,pn+h*kp1/2)
 g=g./norm(g)
 kt2=g[1]
 kp2=g[2]
 g=grad(tn+h*kt2/2,pn+h*kp2/2)
 g=g./norm(g)
 kt3=g[1]
 kp3=g[2]
 g=grad(tn+h*kt3,pn+h*kp3)
 g=g./norm(g)
 kt4=g[1]
 kp4=g[2]


 tt=tn+h/6*(kt1+2*kt2+2*kt3+kt4)
 pp=pn+h/6*(kp1+2*kp2+2*kp3+kp4)

if tt<0
tt=0
pp=pp+pi
end
if tt>pi
tt=pi
pp=pp+pi
end

ncp=length(cps[:,1])
for j=1:ncp
if norm([tt pp]-cps[j,:]')<h
check=false
tt=cps[j,1]
pp=cps[j,2]
break
end
end

gp2=[gp2;tt pp]
tn=tt
pn=pp
end

gp1=gp1[2:end,:]
gp1=reverse(gp1,dims=1)

gp=[gp1;gp2]

if flag==1
out=gp
else
out=0
end

return out
end
#############################################################
#############################################################
function getpoints(n)

ts=collect(range(0,2*pi,n+1))
ts=ts[1:end-1]
r=2
pts=zeros(n,2)
cnt=1

for i=1:n
 pts[i,:]=[pi/2+.5*r*cos(ts[i]) pi+r*sin(ts[i])]
end

pts=[pts; 0 pi; pi pi]

return pts
end
#############################################################
#############################################################
function plotallpaths2d(n)
f(t,p)=cos(p)*sin(t)+cos(2*t)
pts=getpoints(n)

p = range(0, 2*pi, length=100)
t = range(0, pi, length=100)

z = @. f(t', p)


fig=contour(t, p, z, levels=40, color=:plasma, xlabel="theta", ylabel="phi")


gp=gradpath(pts[end-1,1],pts[end-1,2])
plot!(gp[:,1],gp[:,2],legend=false)
scatter!([pts[end-1,1]],[pts[end-1,2]],legend=false)
gp=gradpath(pts[end,1],pts[end,2])
plot!(gp[:,1],gp[:,2],legend=false)
scatter!([pts[end,1]],[pts[end,2]],legend=false)

cnt=0
np=length(pts[:,1])-2
for i=1:np
gp=gradpathflag(pts[i,1],pts[i,2])
if gp==0
cnt=cnt+1
scatter!([pts[i,1]],[pts[i,2]],legend=false)
else
plot!(gp[:,1],gp[:,2],legend=false)
scatter!([pts[i,1]],[pts[i,2]],legend=false)
end
end
display("number of degenerate paths skipped")
display(cnt)

return fig
end
#############################################################
#############################################################
function plotallpaths2drk4(n)
f(t,p)=cos(p)*sin(t)+cos(2*t)
pts=getpoints(n)

p = range(0, 2*pi, length=100)
t = range(0, pi, length=100)

z = @. f(t', p)


fig=contour(t, p, z, levels=40, color=:plasma, xlabel="theta", ylabel="phi")


gp=gradpath(pts[end-1,1],pts[end-1,2])
plot!(gp[:,1],gp[:,2],legend=false)
scatter!([pts[end-1,1]],[pts[end-1,2]],legend=false)
gp=gradpathrk4(pts[end,1],pts[end,2])
plot!(gp[:,1],gp[:,2],legend=false)
scatter!([pts[end,1]],[pts[end,2]],legend=false)

cnt=0
np=length(pts[:,1])-2
for i=1:np
gp=gradpathflagrk4(pts[i,1],pts[i,2])
if gp==0
cnt=cnt+1
scatter!([pts[i,1]],[pts[i,2]],legend=false)
else
plot!(gp[:,1],gp[:,2],legend=false)
scatter!([pts[i,1]],[pts[i,2]],legend=false)
end
end
display("number of degenerate paths skipped")
display(cnt)

return fig
end
#############################################################
#############################################################
function plotallpaths3d(n)

f(t,p)=cos.(p).*sin.(t)+cos.(2 .*t)

pts=getpoints(n)

n = 300
t = LinRange(0,pi,n)
p = LinRange(0,2*pi,n)

tt = t*ones(n)'
pp = ones(n)*p'

xx = @. sin(t)*cos(p')
yy = @. sin(t)*sin(p')
zz = @. cos(t)*ones(n)'

fig=surface(xx, yy, zz, fill_z = f(tt,pp) , size=(600,600))

gp=gradpath(pts[end-1,1],pts[end-1,2])
gp=spheretocart(gp)
plot!(gp[:,1], gp[:,2], gp[:,3], legend=false, linewidth=5)
gp=gradpath(pts[end,1],pts[end,2])
gp=spheretocart(gp)
plot!(gp[:,1], gp[:,2], gp[:,3], legend=false, linewidth=5)


cnt=0
np=length(pts[:,1])-2
for i=1:np
gp=gradpathflag(pts[i,1],pts[i,2])
if gp==0
cnt=cnt+1
else
gp=spheretocart(gp)
plot!(gp[:,1], gp[:,2], gp[:,3], legend=false, linewidth=5)
end
end
display("number of degenerate paths skipped")
display(cnt)


return fig
end
#############################################################
#############################################################
function plotallpaths3drk4(n)

f(t,p)=cos.(p).*sin.(t)+cos.(2 .*t)

pts=getpoints(n)

n = 300
t = LinRange(0,pi,n)
p = LinRange(0,2*pi,n)

tt = t*ones(n)'
pp = ones(n)*p'

xx = @. sin(t)*cos(p')
yy = @. sin(t)*sin(p')
zz = @. cos(t)*ones(n)'

fig=surface(xx, yy, zz, fill_z = f(tt,pp) , size=(600,600))

gp=gradpathrk4(pts[end-1,1],pts[end-1,2])
gp=spheretocart(gp)
plot!(gp[:,1], gp[:,2], gp[:,3], legend=false, linewidth=5)
gp=gradpathrk4(pts[end,1],pts[end,2])
gp=spheretocart(gp)
plot!(gp[:,1], gp[:,2], gp[:,3], legend=false, linewidth=5)


cnt=0
np=length(pts[:,1])-2
for i=1:np
gp=gradpathflagrk4(pts[i,1],pts[i,2])
if gp==0
cnt=cnt+1
else
gp=spheretocart(gp)
plot!(gp[:,1], gp[:,2], gp[:,3], legend=false, linewidth=5)
end
end
display("number of degenerate paths skipped")
display(cnt)


return fig
end
#############################################################
#############################################################
