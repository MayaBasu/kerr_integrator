import kerrgeopy as kg
from math import cos, pi
import numpy as np
import matplotlib.pyplot as plt
#orbit = kg.StableOrbit(0.0, 12, 0.9, 0.5)
orbit = kg.StableOrbit.from_constants(0.99,0.877,1.9,1.2)

print(orbit.constants_of_motion())
t, r, theta, phi = orbit.trajectory()
#print(kg.constants_of_motion(0.0, 12, 0.9, 1.0))

#fig, ax = orbit.plot(0,1)
print(orbit.x)

time = np.linspace(0,5,200)

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




with open('../outpu3t.json', 'r') as file:
    streamdata = json.load(file)["h"]
print(streamdata)

datap = data["phi"]
datat = data["theta"]
datar = data["radial"]
print(datar)
# Print the data
fig = plt.figure(figsize = (6, 6))
ax = plt.axes(projection = '3d')

fig2 = plt.figure()
ax2 = plt.axes()
x = []
yr = []
yt = []
yp = []
z = []
x3 = []
y3 = []
z3 =[]


streamdatalist = []
streamdatax = []
for i in range(len(datar)):
    r = datar[i][1]
    theta = datat[i][1]
    phi = datap[i][1]
    x.append(datat[i][0])
    yr.append(r)
    yp.append(phi)
    yt.append(theta)

    x3.append(r*m.sin(theta)*m.cos(phi))
    y3.append(r*m.sin(theta)*m.sin(phi))
    z3.append(r*m.cos(theta))
    streamdatalist.append(streamdata[i][1])
    streamdatax.append(streamdata[i][0])
    print(r)
ax.plot3D(x3, y3, z3, 'blue')

plt.figure(figsize=(20,4))

plt.subplot(1,4,1)
plt.plot(x, yr)
plt.xlabel("$\lambda$")
plt.ylabel(r"$t(\lambda)$")

plt.subplot(1,4,2)
plt.plot(x ,yt)
plt.xlabel("$\lambda$")
plt.ylabel("$r(\lambda)$")

plt.subplot(1,4,3)
plt.plot(x, yp)
plt.xlabel("$\lambda$")
plt.ylabel(r"$\theta(\lambda)$")

plt.figure()
plt.plot(streamdatax,streamdatalist)




#print(dataarray)
#print(data["data"])




print(kg.apex_from_constants(0.99, 0.99, 6.0, 3))
plt.show()