from mpl_toolkits import mplot3d

import math
import numpy as np
import matplotlib.pyplot as plt
import collections
from collections import Counter

#free space path loss

freq=5180e6 # in Hz
c=3e8 # speed of light in vaccum in m/s
Pt= 20 #power transmitted dBm
SNR= 30 #desired SNR
noise= -85 #noise floor -85 dBm
step=1 #in the points search

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
p1=[min(xymin),ymin,desiredAltitude]
p2=[max(xymin),ymin,desiredAltitude]
p3=[min(xymax),ymax,desiredAltitude]
p4=[max(xymax),ymax,desiredAltitude]

#define lines 
lz=[desiredAltitude,desiredAltitude]
#between p1 and p2
l1x=[p1[0],p2[0]]
l1y=[p1[1],p2[1]]
#between p3 and p4
l2x=[p3[0],p4[0]]
l2y=[p3[1],p3[1]]
#between p1 and ideal
l3x=[p1[0],idealPos[0]]
l3y=[p1[1],idealPos[1]]
#between p2 and ideal
l4x=[p2[0],idealPos[0]]
l4y=[p2[1],idealPos[1]]
#between p3 and ideal
l5x=[p3[0],idealPos[0]]
l5y=[p3[1],idealPos[1]]
#between p4 and ideal
l6x=[p4[0],idealPos[0]]
l6y=[p4[1],idealPos[1]]
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


plt.show()