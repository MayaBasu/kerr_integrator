import kerrgeopy as kg
import numpy as np
import matplotlib.pyplot as plt
import math as m
import json

time = np.linspace(0, 40, 200000)
orbit = kg.StableOrbit.from_constants(0.999, 0.87, 1.9, 1.26)
print("a p e x")
print(kg.apex_from_constants(0.999, 0.87, 1.9, 1.26))
t, r_kg, theta, phi = orbit.trajectory()

print('a   E   Lz    Q')
print(orbit.constants_of_motion())


print(" initial r is ")
print(r_kg(time))

fig, axs = plt.subplots(3, 3)

axs[0, 0].plot(time, r_kg(time))
axs[0, 0].set_xlabel("$\lambda$")
axs[0, 0].set_ylabel("$r(\lambda)$")

axs[0, 1].plot(time, theta(time))
axs[0, 1].set_xlabel("$\lambda$")
axs[0, 1].set_ylabel(r"$\theta(\lambda)$")

axs[0, 2].plot(time, phi(time))
axs[0, 2].set_xlabel("$\lambda$")
axs[0, 2].set_ylabel(r"$\phi(\lambda)$")

with open('../ballistic_graph.json', 'r') as file:
    ballistic_data = json.load(file)

with open('../stream_width.json', 'r') as file:
    stream_data = json.load(file)["h"]


phi_graph = ballistic_data["phi"]
theta_graph = ballistic_data["theta"]
radial_graph = ballistic_data["radial"]

phi_values = []
theta_values = []
radial_values = []

x3 = []
y3 = []
z3 = []


r_ratio = []

x_axis = []
stream_height = []
stream_x_axis = []



three_d_fig = plt.figure(figsize=(6, 6))
ax = plt.axes(projection='3d')


for i in range(len(radial_graph)):
    r = radial_graph[i][1]
    theta = theta_graph[i][1]
    phi = phi_graph[i][1]

    radial_values.append(r)
    theta_values.append(theta)
    phi_values.append(phi)

    x3.append(r * m.sin(theta) * m.cos(phi))
    y3.append(r * m.sin(theta) * m.sin(phi))
    z3.append(r * m.cos(theta))

    stream_height.append(stream_data[i][1])
    stream_x_axis.append(r)

    x_axis.append(radial_graph[i][0])


for i in range(len(x_axis)):

    r_kg_value = r_kg(x_axis[i])
    r_ratio.append(radial_graph[i][1]/(r_kg_value))


print(r_ratio)
ax.plot3D(x3, y3, z3, 'blue')


axs[1, 0].plot(x_axis, radial_values)
plt.xlabel("$\lambda$")
plt.ylabel("$r(\lambda)$")

axs[1, 1].plot(x_axis, theta_values)
plt.xlabel("$\lambda$")
plt.ylabel(r"$\theta(\lambda)$")

axs[1, 2].plot(x_axis, phi_values)
plt.xlabel("$\lambda$")
plt.ylabel(r"$\phi(\lambda)$")

axs[2,1].plot(x_axis,r_ratio)


axs[2,2].plot(x3,y3)
axs[2,2].scatter(x3[200],y3[200])

#plt.plot(stream_x_axis, r_ratio)

fig.tight_layout()
plt.show()
