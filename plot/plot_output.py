import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

data = np.loadtxt('data0', skiprows=0)


x = data[:,0];
y = data[:,1];
z = data[:,2];
#z = np.reshape(data[:,3],(100,100));

ax.plot_trisurf(x, y, z);


plt.show()

