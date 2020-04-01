# UAV Power Simulator
The simulator presented here determines the trajectory for a UAV serving as a relay in a network. Then a comparison between the power consumed in hovering and the power consumed in trajectory is given.

The positioning of the UAVs acting as FAPs in the flying network is given in FAPs.txt. In line 2 the number of FAPs is given, then after "Positions(x,y,z):", the following lines include the coordinates for each of the FAPs. The last line corresponds to the traffic demands (Mbit/s) for each of the FAPs.

To calculate the free space path loss the following variables are necessary:
  freq=5180e6 # in Hz
  c=3e8 # speed of light in vaccum in m/s
  Pt= 20 #power transmitted dBm
  noise= -85 #noise floor -85 dBm
  maxMCS=780 # capacity of the shared wireless medium Mbits/s
  
To calculate the power consumption the following variables are necessary (related to UAV specifications and environment):
  rho=1.225 #air density in kg/m³
  W=20 #UAV weight in Newton
  R=0.4 #rotor radius in meters
  A=0.503 #rotor disc area in m²
  omega=300 #blade angular velocity in radians/second
  Utip=120 #tip speed of the rotor blade(m/s)
  d0=0.6 #fuselage drag ratio
  k=0.1 #incremental correction factor to induced power
  V0=4.03 #mean rotor induced velocity in hover
  delta=0.012 #profile drag coefficient
  s=0.05 #rotor solidity
These values should be updated according to the diserid UAV specifications
