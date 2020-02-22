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
step=1

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

xideal=[j[0] for j in pointsArea]
yideal=[j[1] for j in pointsArea]

idealPos=[sum(xideal)/len(pointsArea),sum(yideal)/len(pointsArea),desiredAltitude]

print("Ideal Position="+str(idealPos))

#define points






plt.show()