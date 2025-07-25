import kerrgeopy as kg
import numpy as np
import matplotlib.pyplot as plt
import math as m
import json

time = np.linspace(0, 4, 2000000)
orbit = kg.StableOrbit.from_constants(0.9, 0.9999, 0.5*6.5, 6.5*6.5*(1.0-0.5*0.5))
print("a p e x")

t_kg, r_kg, theta_kg, phi_kg = orbit.trajectory()

print('a   E   Lz    Q')
print(orbit.constants_of_motion())


print(" initial r is ")
#print(r_kg(time))

#trim r_kg
initial_index = 0
found = False

#for index in range(0,20):
 #   if not found:
  #      if r_kg(4*index/20)>3:
    #        initial_index = index
   #         found = True

print(r_kg(0),theta_kg(0))
#time = np.linspace(4*initial_index/200000, 4, 200000)

print("THe INTDEX IS " + str(initial_index))
print("with value "+str(r_kg(4*initial_index/200000)))
print("the INITIAL THETA IS AT " + str(theta_kg(4 * initial_index / 200000)))
fig, axs = plt.subplots(3, 4)

fig.tight_layout()
axs[0, 0].plot(time, r_kg(time))
axs[0, 0].set_xlabel("$\lambda$")
axs[0, 0].set_ylabel("$r(\lambda)$")

axs[0, 1].plot(time, theta_kg(time))
axs[0, 1].set_xlabel("$\lambda$")
axs[0, 1].set_ylabel(r"$\theta(\lambda)$")

axs[0, 2].plot(time, phi_kg(time))
axs[0, 2].set_xlabel("$\lambda$")
axs[0, 2].set_ylabel(r"$\phi(\lambda)$")

axs[0, 3].plot(time, t_kg(time))
axs[0, 3].set_xlabel("$\lambda$")
axs[0, 3].set_ylabel(r"$t(\lambda)$")

with open('../maingraph.txt', 'r') as file:
    ballistic_data = json.load(file)


phi_graph = ballistic_data["phi_graph"]
theta_graph = ballistic_data["theta_graph"]
radial_graph = ballistic_data["radial_graph"]
t_graph = ballistic_data["t_graph"]
#intersection_points = ballistic_data["possible_self_intersections"]
stream_data = ballistic_data["stream_height"]



phi_values = []
theta_values = []
radial_values = []
t_values = []

x3 = []
y3 = []
z3 = []


r_ratio = []
theta_ratio = []
phi_ratio = []
t_ratio = []

x_axis = []
stream_height = []
stream_x_axis = []








for i in range(len(radial_graph)):
    r = radial_graph[i]
    theta = theta_graph[i]
    phi = phi_graph[i]
    t = t_graph[i]

    radial_values.append(r)
    theta_values.append(theta)
    phi_values.append(phi)
    t_values.append(t)

    x3.append(r * m.sin(theta) * m.cos(phi))
    y3.append(r * m.sin(theta) * m.sin(phi))
    z3.append(r * m.cos(theta))

    stream_height.append(stream_data[i])
    stream_x_axis.append(r)

    x_axis.append(i*0.01)

initial_condition_index = 0
initial_passage = False
second_out = False
found = False


print("THe INTDEX IS " + str(initial_condition_index))
print("with value "+str(radial_values[initial_condition_index]))
print("the INITIAL THETA IS AT " + str(theta_values[initial_condition_index]))
print('the stream at this point is' + str(stream_height[initial_condition_index]))
print("The derivative is " + str((stream_height[initial_condition_index]-stream_height[initial_condition_index-1])/(radial_values[initial_condition_index]-radial_values[initial_condition_index-1])))

initial_index = 0

for i in range(len(x_axis)):

    r_kg_value = r_kg(x_axis[i]+4*initial_index/200000)
    if r_kg_value > 0.00001:
        r_ratio.append(radial_graph[i]/(r_kg_value))
    else:
        print("avoided divide by zero")
        r_ratio.append(1)

    theta_kg_value = theta_kg(x_axis[i] + 4 * initial_index / 200000)
    theta_ratio.append(theta_graph[i]/(theta_kg_value))

    phi_kg_value = phi_kg(x_axis[i] + 4 * initial_index / 200000)
    if phi_kg_value > 0.00001:
        phi_ratio.append(phi_graph[i]/(phi_kg_value))
    else:
        phi_ratio.append(1)
        print("avoided divide by zero")

    t_kg_value = t_kg(x_axis[i])
    if t_kg_value > 0.00001:
        t_ratio.append(t_graph[i]/t_kg_value)
    else:
        print("avoided divide by zero")
        t_ratio.append(1.0)


#print(theta_ratio)



#print(r_ratio)


print(theta_values)
axs[1, 0].plot(x_axis, radial_values)
plt.xlabel("$\lambda$")
plt.ylabel("$r(\lambda)$")

axs[1, 1].plot(x_axis, theta_values)
plt.xlabel("$\lambda$")
plt.ylabel(r"$\theta(\lambda)$")

axs[1, 2].plot(x_axis, phi_values)
plt.xlabel("$\lambda$")
plt.ylabel(r"$\phi(\lambda)$")

axs[1, 3].plot(x_axis, t_values)
plt.xlabel("$\lambda$")
plt.ylabel(r"$t(\lambda)$")

axs[2,0].plot(x_axis,r_ratio)
axs[2,1].plot(x_axis,theta_ratio)
axs[2,2].plot(x_axis,phi_ratio)
axs[2,3].plot(x_axis,t_ratio)

stream_figure = plt.figure()
axxis = plt.axes()
axxis.plot(stream_x_axis,stream_height)
axxis.scatter(stream_x_axis[initial_condition_index],stream_height[initial_condition_index])

fig = plt.figure(figsize=(6, 6))
azis = plt.axes()

azis.plot(x3,y3)
#for i, intersection in enumerate(intersection_points):

 #   color = (0,1,i*(1/(len(intersection_points))))
 #   azis.scatter(x3[intersection[0]],y3[intersection[0]], color = color)
 #   azis.scatter(x3[intersection[0]+1],y3[intersection[0]+1], color = color)

  #  color = (1,0,i*(1/(len(intersection_points))))

  #  azis.scatter(x3[intersection[1]+1],y3[intersection[1]+1], color = color)
  #  azis.scatter(x3[intersection[1]],y3[intersection[1]], color = color)
#


#axs[2,2].scatter(x3[481],y3[481])

#axs[2,2].scatter(x3[382],y3[382])

#plt.plot(stream_x_axis, r_ratio)

three_d_fig = plt.figure(figsize=(6, 6))
ax = plt.axes(projection='3d')
ax.plot3D(x3, y3, z3, 'blue')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

fig.tight_layout()
plt.show()
