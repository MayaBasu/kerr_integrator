import kerrgeopy as kg
from math import cos, pi
import numpy as np
import matplotlib.pyplot as plt
#orbit = kg.StableOrbit(0.0, 12, 0.9, 0.5)
orbit = kg.StableOrbit.from_constants(0.99,0.89,1.9,1.2)
print(orbit.constants_of_motion())
t, r, theta, phi = orbit.trajectory()
#print(kg.constants_of_motion(0.0, 12, 0.9, 1.0))

#fig, ax = orbit.plot(0,1)


time = np.linspace(0,50,200)

plt.figure(figsize=(20,4))

plt.subplot(1,4,1)
plt.plot(time, t(time))
plt.xlabel("$\lambda$")
plt.ylabel(r"$t(\lambda)$")

plt.subplot(1,4,2)
plt.plot(time, r(time))
plt.xlabel("$\lambda$")
plt.ylabel("$r(\lambda)$")

plt.subplot(1,4,3)
plt.plot(time, theta(time))
plt.xlabel("$\lambda$")
plt.ylabel(r"$\theta(\lambda)$")

plt.subplot(1,4,4)
plt.plot(time, phi(time))
plt.xlabel("$\lambda$")
plt.ylabel(r"$\phi(\lambda)$")



import numpy as np
from mpl_toolkits import mplot3d
import math as m
print("hello")

import json

with open('../outpu2t.json', 'r') as file:
    data = json.load(file)
datap = data["phi"]
datat = data["theta"]
datar = data["radial"]
# Print the data
#fig = plt.figure(figsize = (6, 6))
#ax = pyplot.axes(projection = '3d')

fig2 = plt.figure()
ax2 = plt.axes()
x = []
y = []
z = []
for i in range(len(datar)):
    r = datar[i][1]
    theta = datat[i][1]
    phi = datap[i][1]
    x.append(datar[i][0])
    y.append(datar[i][1])

#   x.append(r*m.sin(theta)*m.cos(phi))
#  y.append(r*m.sin(theta)*m.sin(phi))
# z.append(r*m.cos(theta))
#ax.plot3D(x, y, z, 'blue')
ax2.plot(x,y)
print(x)

print(max(datap))


#print(dataarray)
#print(data["data"])




print(kg.apex_from_constants(0.99, 0.99, 6.0, 3))
plt.show()