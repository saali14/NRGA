from plyfile import PlyData, PlyElement
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np


finalPath_Line = '/home/ali/02_PHD_Ali/03_INDEPENDENT_RESEARCH/02_NRGA/NRGA_Application/DATA/INPUT/'
## We have to show Gravitational Potential Energy The Points of the 'Line' dataset!
plydata = PlyData.read(finalPath_Line + 'CosineSurface.ply')
x_line = plydata['vertex']['x']
y_line = plydata['vertex']['y']
z_line = plydata['vertex']['z']
numPtsRef = np.column_stack((x_line, y_line, z_line))


fig = plt.figure()
ax = fig.gca(projection='3d')


def Gravity(x, y, z): 	
	G   = 0.4
    	eps = 0.1
    	dt  = 0.006
	#fx, fy, fz = np.meshgrid(np.arange(-2.0, 2.0, 8), np.arange(-2.0, 2.0, 8), np.arange(-2.0, 2.0, 8))
	fx = [(x - numPtsRef[i][0])/ np.sqrt((x - numPtsRef[i][0])**2 + (y - numPtsRef[i][1])**2 + (z - numPtsRef[i][2])**2 + eps**2) * (G / ((x - numPtsRef[i][0])**2 + (y - numPtsRef[i][1])**2 + (z - numPtsRef[i][2])**2 + eps**2)) * (dt**2) for i in range(len(numPtsRef))]
	fy = [(y - numPtsRef[i][1])/ np.sqrt((x - numPtsRef[i][0])**2 + (y - numPtsRef[i][1])**2 + (z - numPtsRef[i][2])**2 + eps**2) * (G / ((x - numPtsRef[i][0])**2 + (y - numPtsRef[i][1])**2 + (z - numPtsRef[i][2])**2 + eps**2)) * (dt**2) for i in range(len(numPtsRef))]
	fz = [(z - numPtsRef[i][2])/ np.sqrt((x - numPtsRef[i][0])**2 + (y - numPtsRef[i][1])**2 + (z - numPtsRef[i][2])**2 + eps**2) * (G / ((x - numPtsRef[i][0])**2 + (y - numPtsRef[i][1])**2 + (z - numPtsRef[i][2])**2 + eps**2)) * (dt**2) for i in range(len(numPtsRef))]
	
	dispX = x - np.sum(fx) 
	dispY = y - np.sum(fy)
	dispZ = z - np.sum(fz) 

	return dispX, dispY, dispZ 

# prepare a meshgrid 
x, y, z = np.meshgrid(np.arange(-4.0, 4.0, 0.5),
     		      np.arange(-4.0, 4.0, 0.5),
     		      np.arange(-4.0, 4.0, 0.5))


u, v, w = Gravity(x,y,z)

ax.quiver(u, v, w, x, y, z, length=0.05)
ax.set_title('Varying Density')
#ax.quiver(u, v, w, x, y, z, length=0.05)


plt.show()

