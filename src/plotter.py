import kerrgeopy as kg
import numpy as np
import matplotlib.pyplot as plt
import math as m
import json

fig, axs = plt.subplots(3, 4)

with open('../ballistic_graph.json', 'r') as file:
    ballistic_data = json.load(file)


phi_graph = ballistic_data["phi_graph"]
theta_graph = ballistic_data["theta_graph"]
radial_graph = ballistic_data["radial_graph"]
t_graph = ballistic_data["t_graph"]
stream_data = ballistic_data["stream_height"]
x_axis = np.linspace(0,4,4000)



axs[1, 0].plot(x_axis, radial_graph)
plt.xlabel("$\lambda$")
plt.ylabel("$r(\lambda)$")

axs[1, 1].plot(x_axis, theta_graph)
plt.xlabel("$\lambda$")
plt.ylabel(r"$\theta(\lambda)$")

axs[1, 2].plot(x_axis, phi_graph)
plt.xlabel("$\lambda$")
plt.ylabel(r"$\phi(\lambda)$")

axs[1, 3].plot(x_axis, t_graph)
plt.xlabel("$\lambda$")
plt.ylabel(r"$t(\lambda)$")

stream_figure = plt.figure()
axxis = plt.axes()
axxis.plot(radial_graph,stream_data)


fig.tight_layout()
plt.show()
