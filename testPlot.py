from mpl_toolkits import mplot3d

import math
import numpy as np
import matplotlib.pyplot as plt
import collections
from collections import Counter
from scipy import optimize

#free space path loss

freq=5180e6 # in Hz
c=3e8 # speed of light in vaccum in m/s
Pt= 20 #power transmitted dBm
SNR= 30 #desired SNR
noise= -85 #noise floor -85 dBm
step=1 #in the points search

#variables to calculate energy consumption
rho=1.225
W=20
R=0.4
A=0.503
omega=300
Utip=120
d0=0.6
k=0.1
V0=4.03
delta=0.012
s=0.05

exponent= (-SNR-noise+Pt+20*math.log10(c/(4*freq*math.pi)))/20

distance = math.pow(10, exponent) #radius for which the power recieved is equal or greater than the desired

print('radius='+str(distance))


#figure for ranges
fig = plt.figure(1)

ax = plt.axes(projection='3d')
ax.set_title('SNR for FMAPs')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

x= [0,30] #positions of the FMAPs in x
y= [0,10] #positions of the FMAPs in y
z= [0,20] #positions of the FMAPs in z


#generate FMAPs
i=0
while i < len(x):
    ax.scatter(x[i],y[i],z[i],marker='o', label="FMAP"+str(i+1))
    i+=1

#generate spheres
i=0
while i<len(x):
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    xs = x[i] + distance * np.outer(np.cos(u), np.sin(v))
    ys = y[i] + distance* np.outer(np.sin(u), np.sin(v))
    zs = z[i] + distance * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(xs, ys, zs,  rstride=4, cstride=4)
 
    i+=1

leg=ax.legend()


#figure for desired area
fig = plt.figure(2)

ax = plt.axes(projection='3d')
ax.set_title('Intersection of SNR>=threshold')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

#calculate points where all of the SNR is >=20dB

i=0

pd= [[]]*len(x)
dist = [None]*len(x)
print (dist)

while i<len(x):
    pd[i]= np.array((x[i],y[i],z[i]))
    i+=1

print(pd)

xmax,ymax,zmax=max(x),max(y),max(z)

xd,yd,zd=0,0,0

validPoints = []

while xd <= xmax:
    yd=0
    count=0
    while yd <= ymax:
        zd=0
        count=0
        while zd <= zmax:
            currentPoint=np.array((xd,yd,zd))
            print('Current Point ='+str(currentPoint))
            i=0
            count=0
            while i<len(x):
                dist[i] = np.linalg.norm(pd[i]-currentPoint)
                if(dist[i]>0.0):
                    Pr=Pt+20*math.log10(c/(4*freq*dist[i]*math.pi))
                    print(Pr)
                    if((Pr-noise)>=SNR):
                        count+=1
                dist[i]=None
                i+=1    
            if(count==len(x)):
                validPoints.append(currentPoint)
            zd+=step
        yd+=step
    xd+=step

#plot for the points for the volume admissible

validAltitudes=[]

for j in validPoints:
    ax.scatter(j[0],j[1],j[2],marker='o')
    validAltitudes.append(j[2])

occurences= Counter(validAltitudes) #find the altitude with the highest area
desiredAltitude= occurences.most_common(1)[0][0] 

print("Desired Altitude="+str(desiredAltitude))

#plot area for desired altitude

fig = plt.figure(3)

ax = plt.axes(projection='3d')
ax.set_title('Area for Desired Altitude')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

pointsArea=[] #points in the desired area

for j in validPoints:
    if(j[2]==desiredAltitude):
        ax.scatter(j[0],j[1],j[2],marker='o')
        pointsArea.append(j)

#find the centroid for the ideal position

xarray=[j[0] for j in pointsArea]
yarray=[j[1] for j in pointsArea]

idealPos=[sum(xarray)/len(pointsArea),sum(yarray)/len(pointsArea),desiredAltitude]

print("Ideal Position="+str(idealPos))

#define points

ymin=min(yarray)
ymax=max(yarray)

xymin=[] #improve variable names
xymax=[] #same

for j in pointsArea:
    if(j[1]==ymin):
        xymin.append(j[0])
    if(j[1]==ymax):
        xymax.append(j[0])

#trajectory points
point1=[min(xymin),ymin,desiredAltitude]
point2=[max(xymin),ymin,desiredAltitude]
point3=[min(xymax),ymax,desiredAltitude]
point4=[max(xymax),ymax,desiredAltitude]

#define lines 
lz=[desiredAltitude,desiredAltitude]
#between p1 and p2
l1x=[point1[0],point2[0]]
l1y=[point1[1],point2[1]]
#between p3 and p4
l2x=[point3[0],point4[0]]
l2y=[point3[1],point3[1]]
#between p1 and ideal
l3x=[point1[0],idealPos[0]]
l3y=[point1[1],idealPos[1]]
#between p2 and ideal
l4x=[point2[0],idealPos[0]]
l4y=[point2[1],idealPos[1]]
#between p3 and ideal
l5x=[point3[0],idealPos[0]]
l5y=[point3[1],idealPos[1]]
#between p4 and ideal
l6x=[point4[0],idealPos[0]]
l6y=[point4[1],idealPos[1]]
#plot all of them
ax.plot(l1x,l1y,lz)
ax.plot(l2x,l2y,lz)
ax.plot(l3x,l3y,lz)
ax.plot(l4x,l4y,lz)
ax.plot(l5x,l5y,lz)
ax.plot(l6x,l6y,lz)

#trajectory plot

fig = plt.figure(4)

ax = plt.axes(projection='3d')
ax.set_title('Trajectory Line')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.plot(l1x,l1y,lz)
ax.plot(l2x,l2y,lz)
ax.plot(l3x,l3y,lz)
ax.plot(l4x,l4y,lz)
ax.plot(l5x,l5y,lz)
ax.plot(l6x,l6y,lz)

#calculations for energy consumption

P0=(delta/8)*rho*s*A*math.pow(omega,3)*math.pow(R,3)
Pi=(1+k)*(math.pow(W,3/2)/math.sqrt(2*rho*A))

print("P0="+str(P0))
print("P1="+str(Pi))

def P(V):
    firstElement= P0*(1+(3*math.pow(V,2)/(math.pow(Utip,2))))
    square=1+(math.pow(V,4)/(4*math.pow(V0,4)))
    secondElement=Pi*math.pow((math.sqrt(square)-(math.pow(V,2)/(2*math.pow(V0,2)))),1/2)
    thirdElement=(1/2)*d0*rho*s*A*math.pow(V,3)
    return firstElement+secondElement+thirdElement


minVelocity=optimize.fmin(P,0)
minPower=P(minVelocity[0])
powerHover=P(0)


print("Minumum Velocity="+str(minVelocity[0]))
print("Minimum Power="+str(minPower))

#begin calculations for power consumption during trajectory

pn1=np.array(point1)
pn2=np.array(point2)
pn3=np.array(point3)
pn4=np.array(point4)
pnI=np.array(idealPos)


d1=np.linalg.norm(pnI-pn1)
d2=np.linalg.norm(pnI-pn2)
d3=np.linalg.norm(pnI-pn3)
d4=np.linalg.norm(pnI-pn4)
d12=np.linalg.norm(pn2-pn1)
d34=np.linalg.norm(pn4-pn3)


timeTrajectory=(d1/minVelocity)+(d2/minVelocity)+(d3/minVelocity)+(d4/minVelocity)+(d12/minVelocity)+(d34/minVelocity)
powerConsumed=timeTrajectory*minPower+4*powerHover
ifHovering=(timeTrajectory+4)*powerHover

print(timeTrajectory)

print("Trajectory="+str(powerConsumed[0]))
print("Hovering="+str(ifHovering[0]))
print("Compared(trajectory/hovergin):"+str(powerConsumed[0]/ifHovering[0]))



fig = plt.figure(5)
xpart=['Trajectory','Hovering']
ypart=[powerConsumed[0],ifHovering[0]]
y_pos = np.arange(len(ypart))
plt.ylabel('Energy Consumed (Joule)')
plt.bar(y_pos,ypart)
plt.xticks(y_pos,xpart)

#check if longer a trajectory reduces energy consumption compared to hovering

energyV=[]
energyH=[]
xpart=[]
energyRatio=[]
i=1
while i<=30:
    energyV.append(i*minPower+4*powerHover)
    energyH.append((i+4)*powerHover)
    energyRatio.append(energyV[i-1]/energyH[i-1])
    xpart.append(i+4)
    i+=1

fig = plt.figure(6)

plt.plot(xpart,energyV,label="Trajectory")
plt.plot(xpart,energyH,label="Hovering")
plt.xlabel('Time(s)')
plt.ylabel('Energy Consumed(Joule)')
plt.title('Variation of energy consumption')
plt.legend()

fig= plt.figure(7)
plt.plot(xpart,energyRatio)
plt.xlabel('Time(s)')
plt.ylabel('Ratio)')
plt.title('Energy Consumed (Trajectory/Hovering)')


plt.show()