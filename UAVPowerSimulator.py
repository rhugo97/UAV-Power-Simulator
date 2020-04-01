from mpl_toolkits import mplot3d
import sys
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
noise= -85 #noise floor -85 dBm
step=1 #in the points search
maxMCS=780 # capacity of the shared wireless medium Mbits/s

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

f=open("FAPS.txt", "r")

f.readline() #read numbers line
nFMAPS=f.readline()

f.readline() #read positions line

FMAPS=[]
x=[]
y=[]
z=[]
i=0
while i < int(nFMAPS):
    FMAPS.append((f.readline().split(",")))
    x.append(float(FMAPS[i][0]))
    y.append(float(FMAPS[i][1]))
    z.append(float(FMAPS[i][2]))
    i+=1


f.readline() #read traffic line

trafficLine=f.readline().split(",")

i=0
traffic=[]
while i<len(trafficLine):
    traffic.append(float(trafficLine[i]))
    if(traffic[i]>(maxMCS/len(x))):
        traffic[i]=maxMCS/len(x)
    i+=1

def distanceForSNR(SNR):
    exponent= (-SNR-noise+Pt+20*math.log10(c/(4*freq*math.pi)))/20
    return  math.pow(10, exponent) #radius for which the power recieved is equal or greater than the desired


"""
#figure for positions

fig= plt.figure(12)

ax = plt.axes(projection='3d')
ax.set_title('FMAPs Position')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

#generate FAPs
i=0
while i < int(nFMAPS):
    ax.scatter(x[i],y[i],z[i],marker='o', label="FAP"+str(i+1))
    i+=1
"""

#figure for ranges
fig = plt.figure(1)

ax = plt.axes(projection='3d')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

#dictionary for SNR and Data rate relation

dicMCS=[
    {"SNR":11,"data_rate":58.5/len(x)},
    {"SNR":14,"data_rate":117/len(x)},
    {"SNR":18,"data_rate":175.5/len(x)},
    {"SNR":20,"data_rate":234/len(x)},
    {"SNR":24,"data_rate":351/len(x)},
    {"SNR":27,"data_rate":468/len(x)},
    {"SNR":29,"data_rate":526.5/len(x)},
    {"SNR":34,"data_rate":585/len(x)},
    {"SNR":38,"data_rate":702/len(x)},
    {"SNR":40,"data_rate":780/len(x)}
]

#traffic=[200,275] #traffic from the FMAPs in Mbit/s

#map traffic to SNR

data_rate_val=[None]*len(traffic)

j=0
print(len(traffic))
while j < len(traffic):
    i=0
    sums=[]
    while i<10:
        current=dicMCS[i].get("data_rate")-traffic[j]
        if(current>=0):
            sums.append(current)
            if(min(sums)==current):
                data_rate_val[j]=dicMCS[i].get("data_rate")
        i+=1  
    j+=1

SNR_values=[]
for j in data_rate_val:
    i=0
    while i<10:
        data_rate=dicMCS[i].get("data_rate")
        if(j==data_rate):
            SNR_values.append(dicMCS[i].get("SNR"))
        i+=1


print("SNR values= ",SNR_values)
#generate FMAPs
i=0
while i < len(x):
    ax.scatter(x[i],y[i],z[i],marker='o', label="FAP"+str(i+1))
    i+=1

#generate spheres
i=0
while i<len(x):
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    distance=distanceForSNR(SNR_values[i])
    xs = x[i] + distance * np.outer(np.cos(u), np.sin(v))
    ys = y[i] + distance* np.outer(np.sin(u), np.sin(v))
    zs = z[i] + distance * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(xs, ys, zs,  rstride=4, cstride=4,alpha=0.5)
    i+=1

leg=ax.legend()

"""
#figure for desired area
fig = plt.figure(2)

ax = plt.axes(projection='3d')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
"""
#calculate points where all of the SNR is >= threshold

i=0

pd= [[]]*len(x)


while i<len(x):
    pd[i]= np.array((x[i],y[i],z[i]))
    i+=1

print(pd)

#beginning of PREP
xmax,ymax,zmax=max(x),max(y),max(z)
def calculateValidPoints(pd,xmax,ymax,zmax,SNR_values):
    xd,yd,zd=0,0,0
    validPoints = []
    dist = [None]*len(x)
    while xd <= xmax:
        yd=0
        count=0
        while yd <= ymax:
            zd=0
            count=0
            while zd <= zmax:
                currentPoint=np.array((xd,yd,zd))
                #print('Current Point ='+str(currentPoint))
                i=0
                count=0
                while i<len(x):
                    dist[i] = np.linalg.norm(pd[i]-currentPoint)
                    if(dist[i]>0.0):
                        Pr=Pt+20*math.log10(c/(4*freq*dist[i]*math.pi))
                        #print(Pr)
                        if((Pr-noise)>=SNR_values[i]):
                            count+=1
                    #print("FMAP"+str(i)+" with SNR: "+str(SNR_values[i]))        
                    dist[i]=None
                    i+=1    
                if(count==len(x)):
                    validPoints.append(currentPoint)
                zd+=step
            yd+=step
        xd+=step

    return validPoints

happened=0
validPoints= []
while(len(validPoints)==0):
    validPoints=calculateValidPoints(pd,xmax,ymax,zmax,SNR_values)
    print(SNR_values)
    print("Number of Valid Points: "+str(len(validPoints)))
    if(len(validPoints)!=0):
        break
    print("There is no intersection for the current SNR values, therefore there is no valid position for the GW UAV.")
    i=0
    while i<len(SNR_values):
        happened=1
        SNR_values[i]=SNR_values[i]-1
        i+=1

print("SNR values= ",SNR_values)

#plot for the points for the volume admissible

validAltitudes=[]

for j in validPoints:
    #ax.scatter(j[0],j[1],j[2],marker='o')
    validAltitudes.append(j[2])

occurences= Counter(validAltitudes) #find the altitude with the highest area
desiredAltitude= occurences.most_common(1)[0][0] 

print("Desired Altitude="+str(desiredAltitude))

if(happened!=0):
    fig = plt.figure(13)

    ax = plt.axes(projection='3d')
    ax.set_title('SNR for FMAPs')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    i=0
    while i<len(x):
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        distance=distanceForSNR(SNR_values[i])
        xs = x[i] + distance * np.outer(np.cos(u), np.sin(v))
        ys = y[i] + distance* np.outer(np.sin(u), np.sin(v))
        zs = z[i] + distance * np.outer(np.ones(np.size(u)), np.cos(v))
        ax.plot_surface(xs, ys, zs,  rstride=4, cstride=4)
        i+=1

#plot area for desired altitude
"""
fig = plt.figure(3)

ax = plt.axes(projection='3d')
ax.set_title('Area for Desired Altitude')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
"""
pointsArea=[] #points in the desired altitude
for j in validPoints:
    if(j[2]==desiredAltitude):
       # ax.scatter(j[0],j[1],j[2],marker='o')
        pointsArea.append(j)


#find the centroid for the ideal position

xarray=[j[0] for j in pointsArea]
yarray=[j[1] for j in pointsArea]

idealPos=[sum(xarray)/len(pointsArea),sum(yarray)/len(pointsArea),desiredAltitude]

print("Ideal Position="+str(idealPos))

#define points

ymin=min(yarray)
ymax=max(yarray)
xmax=max(xarray)
xmin=min(xarray)

xymin=[] 
xymax=[] 
yxmin=[] 
yxmax=[] 

for j in pointsArea:
    if(j[0]==xmin):
        yxmin.append(j[1])
    if(j[0]==xmax):
        yxmax.append(j[1])
    if(j[1]==ymin):
        xymin.append(j[0])
    if(j[1]==ymax):
        xymax.append(j[0])



def calculateDistanceT(point1,point2,point3,point4,pointI):
    pn1=np.array(point1)
    pn2=np.array(point2)
    pn3=np.array(point3)
    pn4=np.array(point4)
    pnI=np.array(pointI)

    d1=np.linalg.norm(pnI-pn1)
    print("d1=",d1)
    d2=np.linalg.norm(pnI-pn2)
    print("d2=",d2)
    d3=np.linalg.norm(pnI-pn3)
    print("d3=",d3)
    d4=np.linalg.norm(pnI-pn4)
    print("d4=",d4)
    d12=np.linalg.norm(pn2-pn1)
    print("d12=",d12)
    d34=np.linalg.norm(pn4-pn3)
    print("d34=",d34)
    dTotal=d1+d2+d3+d4+d12+d34
    return dTotal

l1x=[0]*3
l1y=[0]*3
l2x=[0]*3
l2y=[0]*3
l3x=[0]*3
l3y=[0]*3
l4x=[0]*3
l4y=[0]*3
l5x=[0]*3
l5y=[0]*3
l6x=[0]*3
l6y=[0]*3

#trajectory 1 points

point1T1=[min(xymin),ymin,desiredAltitude]
point2T1=[max(xymin),ymin,desiredAltitude]
point3T1=[min(xymax),ymax,desiredAltitude]
point4T1=[max(xymax),ymax,desiredAltitude]

#define lines 
lz=[desiredAltitude,desiredAltitude]
#between p1 and p2
l1x[0]=[point1T1[0],point2T1[0]]
l1y[0]=[point1T1[1],point2T1[1]]
#between p3 and p4
l2x[0]=[point3T1[0],point4T1[0]]
l2y[0]=[point3T1[1],point4T1[1]]
#between p1 and ideal
l3x[0]=[point1T1[0],idealPos[0]]
l3y[0]=[point1T1[1],idealPos[1]]
#between p2 and ideal
l4x[0]=[point2T1[0],idealPos[0]]
l4y[0]=[point2T1[1],idealPos[1]]
#between p3 and ideal
l5x[0]=[point3T1[0],idealPos[0]]
l5y[0]=[point3T1[1],idealPos[1]]
#between p4 and ideal
l6x[0]=[point4T1[0],idealPos[0]]
l6y[0]=[point4T1[1],idealPos[1]]
"""#plot all of them
ax.plot(l1x[0],l1y[0],lz)
ax.plot(l2x[0],l2y[0],lz)
ax.plot(l3x[0],l3y[0],lz)
ax.plot(l4x[0],l4y[0],lz)
ax.plot(l5x[0],l5y[0],lz)
ax.plot(l6x[0],l6y[0],lz)
"""
distanceT1=calculateDistanceT(point1T1,point2T1,point3T1,point4T1,idealPos)

#trajectory 2
"""
fig = plt.figure(4)

ax = plt.axes(projection='3d')
ax.set_title('Area for Desired Altitude')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
for j in pointsArea:
    ax.scatter(j[0],j[1],j[2],marker='o')
"""

point1T2=[xmin,min(yxmin),desiredAltitude]
point2T2=[xmin,max(yxmin),desiredAltitude]
point3T2=[xmax,min(yxmax),desiredAltitude]
point4T2=[xmax,max(yxmax),desiredAltitude]

distanceT2=calculateDistanceT(point1T2,point2T2,point3T2,point4T2,idealPos)

#define lines 
lz=[desiredAltitude,desiredAltitude]
#between p1 and p2
l1x[1]=[point1T2[0],point2T2[0]]
l1y[1]=[point1T2[1],point2T2[1]]
#between p3 and p4
l2x[1]=[point3T2[0],point4T2[0]]
l2y[1]=[point3T2[1],point4T2[1]]
#between p1 and ideal
l3x[1]=[point1T2[0],idealPos[0]]
l3y[1]=[point1T2[1],idealPos[1]]
#between p2 and ideal
l4x[1]=[point2T2[0],idealPos[0]]
l4y[1]=[point2T2[1],idealPos[1]]
#between p3 and ideal
l5x[1]=[point3T2[0],idealPos[0]]
l5y[1]=[point3T2[1],idealPos[1]]
#between p4 and ideal
l6x[1]=[point4T2[0],idealPos[0]]
l6y[1]=[point4T2[1],idealPos[1]]
"""#plot all of them
ax.plot(l1x[1],l1y[1],lz)
ax.plot(l2x[1],l2y[1],lz)
ax.plot(l3x[1],l3y[1],lz)
ax.plot(l4x[1],l4y[1],lz)
ax.plot(l5x[1],l5y[1],lz)
ax.plot(l6x[1],l6y[1],lz)
"""

#trajectory 3
"""
fig = plt.figure(5)

ax = plt.axes(projection='3d')
ax.set_title('Area for Desired Altitude')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

for j in pointsArea:
    ax.scatter(j[0],j[1],j[2],marker='o')
"""
idealPosXfloor=math.floor(idealPos[0])
idealPosXceil=math.ceil(idealPos[0])
idealPosYfloor=math.floor(idealPos[1])
idealPosYceil=math.ceil(idealPos[1])
yinIdealXfloor=[]
yinIdealXceil=[]
xinIdealYfloor=[]
xinIdealYceil=[]
for j in pointsArea:
    if(j[0]==idealPosXfloor):
        yinIdealXfloor.append(j[1])
    if(j[0]==idealPosXceil):
        yinIdealXceil.append(j[1])

if(len(yinIdealXfloor)>len(yinIdealXceil)):
    xidealPosR=idealPosXfloor
if(len(yinIdealXfloor)<len(yinIdealXceil)):
    xidealPosR=idealPosXceil
if(len(yinIdealXceil)==len(yinIdealXfloor)):
    xidealPosR=round(idealPos[0])

for j in pointsArea:
    if(j[1]==idealPosYfloor):
        xinIdealYfloor.append(j[0])
    if(j[1]==idealPosYceil):
        xinIdealYceil.append(j[0])

if(len(xinIdealYfloor)>len(xinIdealYceil)):
    yidealPosR=idealPosYfloor
if(len(xinIdealYfloor)>len(xinIdealYceil)):
    yidealPosR=idealPosYceil
if(len(xinIdealYceil)==len(xinIdealYfloor)):
    yidealPosR=round(idealPos[1])

    
yinXideal=[]
xinYideal=[]

for j in pointsArea:
    if(j[0]==xidealPosR):
        yinXideal.append(j[1])
    if(j[1]==yidealPosR):
        xinYideal.append(j[0])

point1T3=[xidealPosR,max(yinXideal)]
point2T3=[min(xinYideal),yidealPosR]
point3T3=[max(xinYideal),yidealPosR]
point4T3=[xidealPosR,min(yinXideal)]

idealPosR=[xidealPosR,yidealPosR]

distanceT3=calculateDistanceT(point1T3,point2T3,point3T3,point4T3,idealPosR)

print("Distance T1: ",distanceT1)
print("Distance T2: ",distanceT2)
print("Distance T3: ",distanceT3)

#define lines 
lz=[desiredAltitude,desiredAltitude]
#between p1 and p2
l1x[2]=[point1T3[0],point2T3[0]]
l1y[2]=[point1T3[1],point2T3[1]]
#between p3 and p4
l2x[2]=[point3T3[0],point4T3[0]]
l2y[2]=[point3T3[1],point4T3[1]]
#between p1 and ideal
l3x[2]=[point1T3[0],idealPos[0]]
l3y[2]=[point1T3[1],idealPos[1]]
#between p2 and ideal
l4x[2]=[point2T3[0],idealPos[0]]
l4y[2]=[point2T3[1],idealPos[1]]
#between p3 and ideal
l5x[2]=[point3T3[0],idealPos[0]]
l5y[2]=[point3T3[1],idealPos[1]]
#between p4 and ideal
l6x[2]=[point4T3[0],idealPos[0]]
l6y[2]=[point4T3[1],idealPos[1]]
"""#plot all of them
ax.plot(l1x[2],l1y[2],lz)
ax.plot(l2x[2],l2y[2],lz)
ax.plot(l3x[2],l3y[2],lz)
ax.plot(l4x[2],l4y[2],lz)
ax.plot(l5x[2],l5y[2],lz)
ax.plot(l6x[2],l6y[2],lz)
"""


#trajectory plot for the chosen path
chosen=0
maxDistance=max(distanceT1,distanceT2,distanceT3)

if(maxDistance==distanceT1):
    chosen=0
if(maxDistance==distanceT2):
    chosen=1
if(maxDistance==distanceT3):
    chosen=2


fig = plt.figure(6)

ax = plt.axes(projection='3d')
ax.set_title('Trajectory Line')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.plot(l3x[chosen],l3y[chosen],lz,label="Pc->P1")
ax.plot(l1x[chosen],l1y[chosen],lz,label="P1->P2")
ax.plot(l4x[chosen],l4y[chosen],lz,label="P2->Pc")
ax.plot(l5x[chosen],l5y[chosen],lz,label="Pc->P3")
ax.plot(l2x[chosen],l2y[chosen],lz,label="P3->P4")
ax.plot(l6x[chosen],l6y[chosen],lz,label="P4->Pc")
ax.legend()

print("Trajectory ",chosen+1, "was selected")
#end of PREP
#calculations for energy consumption

P0=(delta/8)*rho*s*A*math.pow(omega,3)*math.pow(R,3)
Pi=(1+k)*(math.pow(W,3/2)/math.sqrt(2*rho*A))

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


print("Total Distance: ", maxDistance)

timeTrajectory=maxDistance/minVelocity
powerConsumed=timeTrajectory*minPower+4*powerHover
ifHovering=(timeTrajectory+4)*powerHover

print("Time="+str(timeTrajectory[0]))

print("Trajectory="+str(powerConsumed[0]))
print("Hovering="+str(ifHovering[0]))
print("Compared(trajectory/hovering):"+str(powerConsumed[0]/ifHovering[0]))



fig = plt.figure(7)
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

fig = plt.figure(8)

plt.plot(xpart,energyV,label="Trajectory")
plt.plot(xpart,energyH,label="Hovering")
plt.xlabel('Time(s)')
plt.ylabel('Energy Consumed(Joule)')
plt.title('Variation of energy consumption')
plt.legend()

fig= plt.figure(9)
plt.plot(xpart,energyRatio)
plt.xlabel('Time(s)')
plt.ylabel('Ratio')
plt.title('Energy Consumed (Trajectory/Hovering)')

#how it impacts the total lifetime of the UAV, assuming that in hovering it lasts 30 min -- needs to be calculated after knowing the exact UAV

timeHovering=1800 #in seconds
totalCapacity=1800*powerHover #amount of joules that can be spent from the battery
numberTrajectory=totalCapacity/powerConsumed[0] #number of times the trajectory is done
print("Number of trajectories done: "+str(numberTrajectory))
totalTrajectory=numberTrajectory*(timeTrajectory[0]+4) #total time the UAV lasts with the trajectory
print("Total time adopting a trajectory="+str(totalTrajectory))

#show changes (in min)

fig= plt.figure(10)
xpart=['Baseline','PREP']
ypart=[timeHovering/60,totalTrajectory/60]
y_pos = np.arange(len(ypart))
plt.ylabel('Total Operational Time of the UAV (min)')
plt.bar(y_pos,ypart)
plt.xticks(y_pos,xpart)


plt.show()
