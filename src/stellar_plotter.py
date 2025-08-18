import math

import kerrgeopy as kg
import numpy as np
import matplotlib.pyplot as plt
import json
import math as m

from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, FFMpegWriter


from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap

colors = {"orange":'#D1603D',"blue":(0.44,0.44,1),"green":"#386150","brown":"#544343"}








def star_chunk_stream_evolution(path):
    fig, stream_axis = plt.subplots(1, 1)
    with open(path, 'r') as file:
        star = json.load(file)
    star_chunks = star["star_chunks"]
    for star_chunk in star_chunks:
        ballistic_data = star_chunk["geodesic_graph"]

        radial_graph = np.array(ballistic_data["radial_graph"])
        stream_data = np.array(ballistic_data["stream_height"])
        stream_axis.plot(radial_graph, stream_data)
def geodesic_stream_evolution(path):
    fig, [stream_axis,close] = plt.subplots(1, 2)
    with open(path, 'r') as file:
        single_ballistic_data = json.load(file)

    ballistic_data = single_ballistic_data

    radial_graph = np.array(ballistic_data["radial_graph"])
    stream_data = np.array(ballistic_data["stream_height"])
    stream_axis.plot(radial_graph, stream_data,color="orange",label="Debris Stream Width")
    stream_axis.set_xlabel("r",size="x-large")
    stream_axis.set_ylabel("Stream Width",size="x-large")
    stream_axis.set_ylim([0,11.5])
    stream_axis.set_xlim([0,10000])
    stream_axis.grid(axis='both')


    close.plot(radial_graph, stream_data,color="orange",label="Debris Stream Width")
    close.set_xlabel("r",size="x-large")
    close.set_ylabel("Stream Width",size="x-large")
    close.set_title("Close up view")
    close.set_xlim([18,25])
    close.set_ylim([0,0.14])
    plt.grid(axis='both')



    fig.suptitle("Debris Stream Width Evolution as a Function of Radial Distance from the SMBH")
    plt.tight_layout()
   # plt.legend()
    plt.savefig("graphs/stream_height")
def compare_components(a, start_search_divisions,path):
    fig, axs = plt.subplots(2, 4, figsize=(9, 5))
    with open(path, 'r') as file:
        star = json.load(file)
    star_chunks = star["star_chunks"]

    for star_chunk in star_chunks:
        ballistic_data = star_chunk["geodesic_graph"]
        num_steps = len(ballistic_data["radial_graph"])
        step_size = ballistic_data["step_size"]
        stellar_params = ballistic_data["stellar_params"]
        E = stellar_params["e"]
        LZ = stellar_params["lz"]
        C = stellar_params["c"]

        phi_graph = np.array(ballistic_data["phi_graph"])
        theta_graph = np.array(ballistic_data["theta_graph"])
        radial_graph = np.array(ballistic_data["radial_graph"])
        t_graph = np.array(ballistic_data["t_graph"])
        print(theta_graph)
        r_phase, theta_phase = find_initial_phases(a, E, LZ, C, radial_graph[0], theta_graph[0], start_search_divisions)
        orbit = kg.StableOrbit.from_constants(a, E, LZ, C, initial_phases=(0, r_phase, theta_phase, 0))
        t_kg, r_kg, theta_kg, phi_kg = orbit.trajectory()

        x_axis = np.linspace(step_size, num_steps * step_size, num_steps)
        color = (0,0,0,0.3)
        axs[0, 0].plot(x_axis, r_kg(x_axis),color = color)
        axs[0, 0].set_xlabel("$\lambda$")
        axs[0, 0].set_ylabel("$r(\lambda)$")

        axs[0, 1].plot(x_axis, theta_kg(x_axis),color = color)
        axs[0, 1].set_xlabel("$\lambda$")
        axs[0, 1].set_ylabel(r"$\theta(\lambda)$")

        axs[0, 2].plot(x_axis, phi_kg(x_axis), color = color)
        axs[0, 2].set_xlabel("$\lambda$")
        axs[0, 2].set_ylabel(r"$\phi(\lambda)$")

        axs[0, 3].plot(x_axis, t_kg(x_axis),color = color)
        axs[0, 3].set_xlabel("$\lambda$")
        axs[0, 3].set_ylabel(r"$t(\lambda)$")


        color = (1,0,1)
        axs[0, 0].plot(x_axis, radial_graph,color = color)
        axs[0, 0].set_xlabel("$\lambda$")
        axs[0, 0].set_ylabel("$r(\lambda)$")

        axs[0, 1].plot(x_axis, theta_graph,color = color)
        axs[0, 1].set_xlabel("$\lambda$")
        axs[0, 1].set_ylabel(r"$\theta(\lambda)$")

        axs[0, 2].plot(x_axis, phi_graph,color = color)
        axs[0, 2].set_xlabel("$\lambda$")
        axs[0, 2].set_ylabel(r"$\phi(\lambda)$")

        axs[0, 3].plot(x_axis, t_graph,color = color)
        axs[0, 3].set_xlabel("$\lambda$")
        axs[0, 3].set_ylabel(r"$t(\lambda)$")





        axs[1, 0].plot(x_axis, radial_graph/r_kg(x_axis))
        axs[1, 0].set_xlabel("$\lambda$")
        axs[1, 0].set_ylabel("$r(\lambda)$")

        axs[1, 1].plot(x_axis, theta_graph/theta_kg(x_axis))
        axs[1, 1].set_xlabel("$\lambda$")
        axs[1, 1].set_ylabel(r"$\theta(\lambda)$")

        axs[1, 2].plot(x_axis, phi_graph/phi_kg(x_axis))
        axs[1, 2].set_xlabel("$\lambda$")
        axs[1, 2].set_ylabel(r"$\phi(\lambda)$")

        axs[1, 3].plot(x_axis, t_graph/t_kg(x_axis))
        axs[1, 3].set_xlabel("$\lambda$")
        axs[1, 3].set_ylabel(r"$t(\lambda)$")
        fig.tight_layout()
        fig.canvas.manager.set_window_title('component_compare')
        fig.suptitle("Comparison of Geodesic code and KerrGeoPy")
        plt.savefig("graphs/component_comparison")
def compare_single_geodesic_components(a, start_search_divisions,path):
    fig, axs = plt.subplots(2, 4, figsize=(9, 5))
    with open(path, 'r') as file:
        single_ballistic_data = json.load(file)


    ballistic_data = single_ballistic_data
    num_steps = len(ballistic_data["radial_graph"])
    step_size = ballistic_data["step_size"]
    stellar_params = ballistic_data["stellar_params"]
    E = stellar_params["e"]
    LZ = stellar_params["lz"]
    C = stellar_params["c"]

    phi_graph = np.array(ballistic_data["phi_graph"])
    theta_graph = np.array(ballistic_data["theta_graph"])
    radial_graph = np.array(ballistic_data["radial_graph"])
    t_graph = np.array(ballistic_data["t_graph"])
    print(theta_graph)
    r_phase, theta_phase = find_initial_phases(a, E, LZ, C, radial_graph[0], theta_graph[0], start_search_divisions)
    orbit = kg.StableOrbit.from_constants(a, E, LZ, C, initial_phases=(0, r_phase, theta_phase, 0))
    t_kg, r_kg, theta_kg, phi_kg = orbit.trajectory()

    x_axis = np.linspace(step_size, num_steps * step_size, num_steps)
    color = (0,0,0,0.3)
    axs[0, 0].plot(x_axis, r_kg(x_axis),color = color,label="KerrGeoPy Baseline")
    axs[0, 0].set_xlabel("$\lambda$")
    axs[0, 0].set_ylabel("$r(\lambda)$")

    axs[0, 1].plot(x_axis, theta_kg(x_axis),color = color)
    axs[0, 1].set_xlabel("$\lambda$")
    axs[0, 1].set_ylabel(r"$\theta(\lambda)$")

    axs[0, 2].plot(x_axis, phi_kg(x_axis), color = color)
    axs[0, 2].set_xlabel("$\lambda$")
    axs[0, 2].set_ylabel(r"$\phi(\lambda)$")

    axs[0, 3].plot(x_axis, t_kg(x_axis),color = color)
    axs[0, 3].set_xlabel("$\lambda$")
    axs[0, 3].set_ylabel(r"$t(\lambda)$")


    color = colors["blue"]
    axs[0, 0].plot(x_axis, radial_graph,color = color,linestyle="dotted",label="Rust Geodesic Integrator")

    axs[0, 0].set_xlabel("$\lambda$")
    axs[0, 0].set_ylabel("$r(\lambda)$")

    axs[0, 1].plot(x_axis, theta_graph,color = color,linestyle="dotted")
    axs[0, 1].set_xlabel("$\lambda$")
    axs[0, 1].set_ylabel(r"$\theta(\lambda)$")

    axs[0, 2].plot(x_axis, phi_graph,color = color,linestyle="dotted")
    axs[0, 2].set_xlabel("$\lambda$")
    axs[0, 2].set_ylabel(r"$\phi(\lambda)$")

    axs[0, 3].plot(x_axis, t_graph,color = color,linestyle="dotted")
    axs[0, 3].set_xlabel("$\lambda$")
    axs[0, 3].set_ylabel(r"$t(\lambda)$")


    for i in range(0,4):

        axs[1, i].plot([x_axis[0],x_axis[-1]],[1,1],color = (0,0,0,0.3))
        axs[1, i].set_ylim([0.99, 1.01])

    color = colors["orange"]
    axs[1, 0].plot(x_axis, radial_graph/r_kg(x_axis),color = color,linestyle="dotted",label="Ratio of outputs")
    fig.legend()
    axs[1, 0].set_xlabel("$\lambda$")
    axs[1, 0].set_ylabel("$r(\lambda)$")
    axs[1, 1].plot(x_axis, theta_graph/theta_kg(x_axis),color = color,linestyle="dotted")
    axs[1, 1].set_xlabel("$\lambda$")
    axs[1, 1].set_ylabel(r"$\theta(\lambda)$")
    axs[1, 2].plot(x_axis, phi_graph/phi_kg(x_axis),color = color,linestyle="dotted")
    axs[1, 2].set_xlabel("$\lambda$")
    axs[1, 2].set_ylabel(r"$\phi(\lambda)$")
    axs[1, 3].plot(x_axis, t_graph/t_kg(x_axis),color = color,linestyle="dotted")
    axs[1, 3].set_xlabel("$\lambda$")
    axs[1, 3].set_ylabel(r"$t(\lambda)$")

    fig.tight_layout()
    fig.canvas.manager.set_window_title('compare_single_geodesic_components')
    fig.suptitle("Comparison of Rust geodesic code and KerrGeoPy")
    plt.savefig("graphs/compare_single_geodesic_components")
def plot_components(path):
    fig, axs = plt.subplots(1, 4, figsize=(9, 5))
    with open(path, 'r') as file:
        star = json.load(file)
    star_chunks = star["star_chunks"]

    for star_chunk in star_chunks:
        ballistic_data = star_chunk["geodesic_graph"]


        phi_graph = np.array(ballistic_data["phi_graph"])
        theta_graph = np.array(ballistic_data["theta_graph"])
        radial_graph = np.array(ballistic_data["radial_graph"])
        t_graph = np.array(ballistic_data["t_graph"])


        axs[0].plot(t_graph, radial_graph)
        axs[0].set_xlabel("$\lambda$")
        axs[0].set_ylabel("$r(\lambda)$")

        axs[ 1].plot(t_graph, theta_graph)
        axs[1].set_xlabel("$\lambda$")
        axs[1].set_ylabel(r"$\theta(\lambda)$")

        axs[2].plot(t_graph, phi_graph)
        axs[2].set_xlabel("$\lambda$")
        axs[2].set_ylabel(r"$\phi(\lambda)$")

        axs[3].plot(t_graph, t_graph)
        axs[3].set_xlabel("$\lambda$")
        axs[3].set_ylabel(r"$t(\lambda)$")

        fig.tight_layout()
        fig.canvas.manager.set_window_title('component_compare')
        fig.suptitle("Comparison of Geodesic code and KerrGeoPy")
        plt.savefig("graphs/component_comparison")
def plot_geodesic_components(mino,path):
    with open(path, 'r') as file:
        single_ballistic_data = json.load(file)
    fig, axs = plt.subplots(1, 4, figsize=(9, 5))
    num_steps = len(single_ballistic_data["radial_graph"])
    step_size = single_ballistic_data["step_size"]

    phi_graph = np.array(single_ballistic_data["phi_graph"])
    theta_graph = np.array(single_ballistic_data["theta_graph"])
    radial_graph = np.array(single_ballistic_data["radial_graph"])
    t_graph = np.array(single_ballistic_data["t_graph"])

    if mino:
        x_axis = np.linspace(step_size, num_steps * step_size, num_steps)
    else:
        x_axis = t_graph

    axs[0].plot(x_axis, radial_graph)
    axs[0].set_xlabel("$\lambda$")
    axs[0].set_ylabel("$r(\lambda)$")

    axs[ 1].plot(x_axis, theta_graph)
    axs[1].set_xlabel("$\lambda$")
    axs[1].set_ylabel(r"$\theta(\lambda)$")

    axs[2].plot(x_axis, phi_graph)
    axs[2].set_xlabel("$\lambda$")
    axs[2].set_ylabel(r"$\phi(\lambda)$")

    axs[3].plot(x_axis, t_graph)
    axs[3].set_xlabel("$\lambda$")
    axs[3].set_ylabel(r"$t(\lambda)$")

    fig.tight_layout()
    fig.canvas.manager.set_window_title('component_compare')
    fig.suptitle("Comparison of Geodesic code and KerrGeoPy")
    plt.savefig("graphs/component_comparison")
def plot_3d_star_chunks(path):
    with open(path, 'r') as file:
        star = json.load(file)
    star_chunks = star["star_chunks"]
    three_d_fig = plt.figure(figsize=(6, 6))
    ax = plt.axes(projection='3d')

    x_ends = []
    y_ends = []
    z_ends = []
    color_list = []
    colors = []


    for (index,star_chunk) in enumerate(star_chunks):

        if index % 10 ==0:

            ballistic_data = star_chunk["geodesic_graph"]


            phi_graph = np.array(ballistic_data["phi_graph"])
            theta_graph = np.array(ballistic_data["theta_graph"])
            radial_graph = np.array(ballistic_data["radial_graph"])

            x = (radial_graph * np.sin(theta_graph) * np.cos(phi_graph))
            y = (radial_graph * np.sin(theta_graph) * np.sin(phi_graph))
            z = (radial_graph * np.cos(theta_graph))
            ax.plot3D(x, y, z,color = (index/len(star_chunks),0,1,0.5))
            x_ends.append(x[-1])
            y_ends.append(y[-1])
            z_ends.append(z[-1])
            #colors.append((index/len(star_chunks),0,1))
            colors.append((index/len(star_chunks),0,1))
            color_list.append(star_chunk["binding_energy"])
 # Red, Green, Blue
    n_bins = 100  # Number of color gradations
    custom_cmap = LinearSegmentedColormap.from_list("custom", colors, N=n_bins)
    scatter = ax.scatter(x_ends,y_ends,z_ends,c=color_list,cmap=custom_cmap)
    cbar = plt.colorbar(scatter,fraction=0.026, pad=0.08)

    cbar.set_label("Specific Binding Energy",size="x-large")

   # ax.set_xlabel("x",size="xx-large")
    ax.set_ylabel("y",size="xx-large")
    ax.set_zlabel("z",size="xx-large")
    plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')
    plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')
    ax.set_title("Angular Debris Spread due to Varying Binding Energies")

    plt.tight_layout()
    plt.savefig("graphs/three_d_trajectories")
def plot_3d_geodesic_animatted(path):
    with open(path, 'r') as file:
        single_ballistic_data = json.load(file)







    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    phi_graph = np.array(single_ballistic_data["phi_graph"])
    theta_graph = np.array(single_ballistic_data["theta_graph"])
    radial_graph = np.array(single_ballistic_data["radial_graph"])
    wrapped_data = single_ballistic_data["wrapped_data"]

    x = (radial_graph * np.sin(theta_graph) * np.cos(phi_graph))
    y = (radial_graph * np.sin(theta_graph) * np.sin(phi_graph))
    z = (radial_graph * np.cos(theta_graph))



    ax.scatter(0,0,0,color="black",label="SMBH location")

    ax.set_xlabel("x",size="xx-large")
    ax.set_ylabel("y",size="xx-large")
    ax.set_zlabel("z",size="xx-large")





    def update(num, x, y, z, line):
        num = num*10000
        line.set_data(x[:num], y[:num])
        line.set_3d_properties(z[:num])
        return line,
    print(len(x))






    line, = ax.plot(x,y,z,color=colors["orange"],label="Center-of-mass geodesic")
    ax.set_title("Elucidated 3-dimensional Orbital Trajectory",size="large")
    plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')

    ani = FuncAnimation(fig, update, frames=50, fargs=(x, y, z, line), interval=10)



    plt.show()





    plt.legend(loc="upper right")




    plt.savefig("graphs/single_three_d_trajectory")
    fig, ax2 = plt.subplots(1, 1, figsize=(3, 3))


    ax2.plot(y,z,color=colors["orange"],label="Projected center-of-mass geodesic")
    ax2.scatter(0,0,color="black",label="SMBH location")

    ax2.set_xlim([-0.2,0.2])
    ax2.set_ylim([-12,12])
    ax2.set_ylabel("z",size="xx-large")
    ax2.set_xlabel("x",size="xx-large")
   # plt.legend()
    plt.tight_layout()
   # ax2.set_title("Close-up Side View of the Stream Precession")


    plt.savefig("graphs/single_three_d_trajectory_projection")



def plot_3d_geodesic(path):
    with open(path, 'r') as file:
        single_ballistic_data = json.load(file)







    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    phi_graph = np.array(single_ballistic_data["phi_graph"])
    theta_graph = np.array(single_ballistic_data["theta_graph"])
    radial_graph = np.array(single_ballistic_data["radial_graph"])
    wrapped_data = single_ballistic_data["wrapped_data"]

    x = (radial_graph * np.sin(theta_graph) * np.cos(phi_graph))
    y = (radial_graph * np.sin(theta_graph) * np.sin(phi_graph))
    z = (radial_graph * np.cos(theta_graph))



    ax.scatter(0,0,0,color="black",label="SMBH location")

    ax.set_xlabel("x",size="xx-large")
    ax.set_ylabel("y",size="xx-large")
    ax.set_zlabel("z",size="xx-large")





    ax.plot(x,y,z,color=colors["orange"],label="Center-of-mass geodesic")
    ax.set_title("Elucidated 3-dimensional Orbital Trajectory",size="large")
    plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')

    plt.legend(loc="upper right")




    plt.savefig("graphs/single_three_d_trajectory")
    fig, ax2 = plt.subplots(1, 1, figsize=(3, 3))


    ax2.plot(y,z,color=colors["orange"],label="Projected center-of-mass geodesic")
    ax2.scatter(0,0,color="black",label="SMBH location")

    ax2.set_xlim([-0.2,0.2])
    ax2.set_ylim([-12,12])
    ax2.set_ylabel("z",size="xx-large")
    ax2.set_xlabel("x",size="xx-large")
    # plt.legend()
    plt.tight_layout()
    # ax2.set_title("Close-up Side View of the Stream Precession")


    plt.savefig("graphs/single_three_d_trajectory_projection")
def plot_stellar_profile(path):
    fraction_masses = []
    binding_energies = []
    with open(path, 'r') as file:
        star = json.load(file)
    star_chunks = star["star_chunks"]

    for star_chunk in star_chunks:
        fraction_of_star = star_chunk["fraction_of_star"]
        binding_energy = star_chunk["binding_energy"]
        fraction_masses.append(fraction_of_star * len(star_chunks))
        binding_energies.append(float(binding_energy))

    stellar_profile, (stellar_profile_axis) = plt.subplots(1, 1, figsize=(6.4, 4), layout="constrained")
    stellar_profile_axis.xaxis.set_inverted(True)
    stellar_profile_axis.plot(binding_energies, fraction_masses)
    stellar_profile.suptitle("Normalized Fraction of Stellar Mass at Each Binding Energy")
    stellar_profile_axis.set_xlim([0.9999,0.9996894])
    stellar_profile_axis.set_xlabel('Specific Binding Energy of Slices')
    plt.xticks(np.linspace(binding_energies[0],binding_energies[-1],5))
    stellar_profile_axis.set_ylabel('Normalized Fraction of Stellar Mass')
    stellar_profile_axis.grid(True)
    stellar_profile.canvas.manager.set_window_title('fraction_with_binding_energy')
    stellar_profile.tight_layout()
    plt.savefig("graphs/fraction_with_binding_energy.png")
def find_initial_phases(a, E, LZ, C, r0, theta0, start_search_divisions):
    orbit = kg.StableOrbit.from_constants(a, E, LZ, C)
    a, p, e, x = kg.apex_from_constants(a, E, LZ, C)
    t_kg, r_kg, theta_kg, phi_kg = orbit.trajectory()

    half_r_period = math.pi / kg.r_frequency(a, p, e, x)
    half_theta_period = math.pi / kg.theta_frequency(a, p, e, x)
    print(theta0)

    initial_r_index = np.argmin(np.absolute(np.array(r_kg(np.linspace(0, half_r_period, start_search_divisions))) - r0))
    initial_theta_index = np.argmin(
        np.absolute(np.array(theta_kg(np.linspace(0, half_theta_period, start_search_divisions))) - theta0))

    r_delay_time = initial_r_index * (half_r_period / start_search_divisions)
    theta_delay_time = initial_theta_index * (half_theta_period / start_search_divisions)

    print("r: Matching the initial condition of {} by taking {} at a delay_time of {}, index is {}".format(r0, r_kg(
        r_delay_time), r_delay_time, initial_r_index))
    print("theta: Matching the initial condition of {} by taking {} at a delay_time of {}, index is {}".format(theta0,
                                                                                                               theta_kg(
                                                                                                                   theta_delay_time),
                                                                                                               theta_delay_time,
                                                                                                               initial_theta_index))

    r_phase = (r_delay_time / half_r_period) * math.pi
    theta_phase = (theta_delay_time / half_theta_period) * math.pi

    return r_phase, theta_phase
def plot_angular_momentum_over_time(path):
    with open(path, 'r') as file:
        star = json.load(file)
    star_chunks = star["star_chunks"]
    fig, axs = plt.subplots(1, 2, figsize=(9, 5))
    angular_momentums = star["angular_momentum_over_time"]
    times = []
    zcomps = []
    xycomps = []
    angles  =[]

    for angular_momentum in angular_momentums:
        zcomp = angular_momentum[3]
        xycomp = m.sqrt(angular_momentum[1]**2+angular_momentum[2]**2)
        total = m.sqrt(zcomp**2+xycomp**2)
        angle = zcomp/total

        angles.append(angle)
        zcomps.append(zcomp)
        xycomps.append(xycomp)


        times.append(angular_momentum[0])



    axs[0].plot(times,zcomps,color=colors["blue"],label="z-component of total angular momentum")
    axs[0].set_xlabel("Time",size="xx-large")
    axs[0].set_ylabel("Specific Angular Momentum",size="x-large")
    axs[0].plot(times,xycomps,color=colors["green"], label="xy-component of total angular momentum")
    axs[1].plot(times[0],angles[0],color=colors["orange"],label="Cosine of approximate total debris inclination")
    axs[1].set_xlabel("Time",size="xx-large")
    axs[1].set_ylabel("$L_z/L$",size="xx-large")

    fig.legend()
    fig.suptitle("Time Evolution of Total Angular Momentum before Initial Self Intersection",size="x-large")
    plt.tight_layout()
def plot_phi_bin(path):
    fig, ax = plt.subplots()


    with open(path, 'r') as file:
        test = json.load(file)

    wrapped_data = test["wrapped_data"]
    for i in [1191]:

        wrapped_data_bin = wrapped_data[i]
        print(wrapped_data_bin)
        xlist =[]
        ylist = []
        colors = []
        color_palatte = ['red','blue','yellow','green','magenta','orange','black','teal','purple','red','blue','yellow','green','magenta','orange','black','teal','purple']

        for wrapped_data_point in wrapped_data_bin:

            r = wrapped_data_point["r"]
            theta = wrapped_data_point["theta"]
            wrap = wrapped_data_point["wrap"]
            phi = wrapped_data_point["phi"]
            h = wrapped_data_point["h"]
            print(h)
            x = r * np.sin(theta) * np.cos(phi)
            y = (r * np.sin(theta) * np.sin(phi))
            xlist.append(x)
            ylist.append(y)
            print("wrap is" + str(wrap))
            colors.append(color_palatte[wrap])

            circle = plt.Circle((x,y ),h , color=(1,0,0,0.1))
            ax.add_patch(circle)

        ax.scatter(xlist,ylist,color=colors)



def plot_end_results_6():

    X = np.array([[(5645140.770046568, -1.524136723443274, -0.8249823924045651, 1.625),(5643336.316784389, -1.5959164991335604, -1.258997188451125, 3.25),(5641551.576797498, -1.449849852922865, -1.1476945004746097, 4.875000000000001)],
                  [(4998456.477494477, 0.2302241507118562, 0.028667614028062362, 1.1180339887498947),(4999666.103926818, 0.05442972490282417, -0.46895129150270687, 2.2360679774997894),(5031907.250613445, -0.10159867001392207, -0.43897580006537384, 3.3541019662496865)],
                  [(5220752.59396204, 0.16060268033130481, 0.08199964528772301, 0.9682458365518551),(5222019.416258693, 0.4106449751683912, 0.39538361666501715, 1.9364916731037103),(5223501.101217638, 0.12979091184032768, -0.08471489192852871, 2.904737509655565)]])
    Y = np.zeros([3,3])
    for i in range(0,len(X)):
        for j in range(0,len(X[0])):
            lx = X[i][j][1]
            ly = X[i][j][2]
            lz = X[i][j][3]
            ltotal = np.sqrt((lx**2+ly**2+lz**2))
            print(lz/ltotal)
            print(X[i][j])

            Y[i][j] = lz/ltotal
    print(Y)
    fig = plt.figure(figsize=(6,6))
    plt.imshow(Y, cmap="inferno",vmin=0.5)
    y_label_list = ['6.2', '4.47', '3.87']
    x_label_list = ['0.25','0.5', '0.75']
    plt.gca().set_xticklabels(x_label_list)
    plt.gca().set_yticklabels(y_label_list)
    plt.gca().set_xticks([0,1,2])
    plt.gca().set_yticks([0,1,2])
    plt.gca().set_xlabel("Initial value of $L_z/L$",size='x-large')
    plt.gca().set_ylabel("Initial Total Angular Momentum",size='x-large')
    #plt.gca().set_xticks([6.5,np.sqrt(20),np.sqrt(15)])

    cb = plt.colorbar(label="Final Total Angular Momentum")
    cb.set_label(label="Final Total Angular Momentum", size='x-large')

    plt.title("Final Angular Momentum Calculation: $10^6 M_{\star}$ SMBH ")
    plt.tight_layout()

def plot_end_results_7():

    X = np.array([[(5643796.555300895, -1.5280668760439735, -0.8265910065507388, 1.625),(5643689.126597131, -1.5957571861942081, -1.253009032464471, 3.25),(5643604.911306672, -1.4414411538593292, -1.1402915921291128, 4.875000000000001)],
                  [(5170670.545151273, 0.2772748202781714, 0.07311343903988908, 1.1180339887498947),(5173111.519487047, 0.0843139215665937, -0.4738239910093439, 2.2360679774997894),(5206750.695044351, -0.08187664850582783, -0.3931909304278861, 3.3541019662496865)],
                  [(5228659.9726113565, 0.096978811028594, 0.20580412998331432, 0.9682458365518551),(5230006.135827178, 0.4521115781950686, 0.2785156952081441, 1.9364916731037103),(5232033.01183234, 0.1505161444651099, -0.0976545607657328, 2.904737509655565)]])


    Y = np.zeros([3,3])
    for i in range(0,len(X)):
        for j in range(0,len(X[0])):
            lx = X[i][j][1]
            ly = X[i][j][2]
            lz = X[i][j][3]
            ltotal = np.sqrt((lx**2+ly**2+lz**2))
            print(lz/ltotal)
            print(X[i][j])

            Y[i][j] = lz/ltotal
    print(Y)
    fig = plt.figure(figsize=(6,6))
    plt.imshow(Y, cmap="inferno",vmin=0.5)
    y_label_list = ['6.2', '4.47', '3.87']
    x_label_list = ['0.25','0.5', '0.75']
    plt.gca().set_xticklabels(x_label_list)
    plt.gca().set_yticklabels(y_label_list)
    plt.gca().set_xticks([0,1,2])
    plt.gca().set_yticks([0,1,2])
    plt.gca().set_xlabel("Initial value of $L_z/L$",size='x-large')
    plt.gca().set_ylabel("Initial Total Angular Momentum",size='x-large')
    #plt.gca().set_xticks([6.5,np.sqrt(20),np.sqrt(15)])

    cb = plt.colorbar()
    cb.set_label(label="Final Total Angular Momentum", size='x-large')

    plt.title("Final Angular Momentum Calculation: $10^7 M_{\star}$ SMBH ")
    plt.tight_layout()


def plot_end_results_8():

    X = np.array([[(5595867.079394745, -1.6246878106488702, -0.8559396778138414, 1.625),(5595895.125356542, -1.5029770284912949, -1.1374841062961698, 4.875000000000001),(5595867.036628354, -1.6694268529187237, -1.267748908627183, 3.25)],
                  [(5650245.741652196, 0.23343981829762459, 0.1384490003171102, 1.1180339887498947),(5652329.939440656, 0.26967128088082326, -0.37705120165216954, 2.2360679774997894),(5685712.3812498385, -0.007118111417913941, -0.3257281742989255, 3.3541019662496865)],
                  [(5652870.777982576, 0.11242753768162313, -0.0584006438888951, 0.9682458365518551),(5653415.193080202, 0.29676004152841173, 0.4427273651136203, 1.9364916731037103),(5654941.812916952, 0.13546579187645824, -0.03824896373817541, 2.904737509655565)]])
    Y = np.zeros([3,3])
    for i in range(0,len(X)):
        for j in range(0,len(X[0])):
            lx = X[i][j][1]
            ly = X[i][j][2]
            lz = X[i][j][3]
            ltotal = np.sqrt((lx**2+ly**2+lz**2))
            print(lz/ltotal)
            print(X[i][j])

            Y[i][j] = lz/ltotal
    print(Y)
    fig = plt.figure(figsize=(6,6))

    plt.imshow(Y, cmap="inferno",vmin=0.5)

    y_label_list = ['6.2', '4.47', '3.87']
    x_label_list = ['0.25','0.5', '0.75']
    plt.gca().set_xticklabels(x_label_list)
    plt.gca().set_yticklabels(y_label_list)
    plt.gca().set_xticks([0,1,2])
    plt.gca().set_yticks([0,1,2])
    plt.gca().set_xlabel("Initial value of $L_z/L$",size='x-large')
    plt.gca().set_ylabel("Initial Total Angular Momentum",size='x-large')
    #plt.gca().set_xticks([6.5,np.sqrt(20),np.sqrt(15)])

    cb = plt.colorbar()
    cb.set_label(label="Final Value of $L_z/L$", size='x-large')

    plt.title("Final Angular Momentum Calculation: $10^8 M_{\star}$ SMBH ")
    plt.tight_layout()





path = "../tiny.txt"

#plot_stellar_profile("stars/star_test1_0.5_6.5_6.txt")


#plot_3d_star_chunks(path)

#plot_end_results_8()
#plot_3d_geodesic_animatted("../test.txt")
plot_geodesic_components(True,path)

plt.show()
