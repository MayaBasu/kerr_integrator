import matplotlib.pyplot as pyplot
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
fig = pyplot.figure(figsize = (6, 6))
#ax = pyplot.axes(projection = '3d')

fig2 = pyplot.figure()
ax2 = pyplot.axes()
x = []
y = []
z = []
for i in range(len(datar)):
    r = datar[i][1]
    theta = datat[i][1]
    phi = datap[i][1]

    x.append(r*m.sin(theta)*m.cos(phi))
    y.append(r*m.sin(theta)*m.sin(phi))
    z.append(r*m.cos(theta))
#ax.plot3D(x, y, z, 'blue')
ax2.plot(datat)
print(max(datap))


#print(dataarray)
#print(data["data"])

pyplot.show()
