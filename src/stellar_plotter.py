import math

import kerrgeopy as kg
import numpy as np
import matplotlib.pyplot as plt
import json
import math as m

from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap

colors = {"orange":'#D1603D',"blue":(0.44,0.44,1),"green":"#386150","brown":"#544343"}

with open('stars/star_0.5_6.5_6.txt', 'r') as file:

    star = json.load(file)
    star_chunks = star["star_chunks"]

with open('maingraph.txt', 'r') as file:
    single_ballistic_data = json.load(file)
def star_chunk_stream_evolution():
    fig, stream_axis = plt.subplots(1, 1)
    for star_chunk in star_chunks:
        ballistic_data = star_chunk["geodesic_graph"]

        radial_graph = np.array(ballistic_data["radial_graph"])
        stream_data = np.array(ballistic_data["stream_height"])
        stream_axis.plot(radial_graph, stream_data)

def geodesic_stream_evolution():
    fig, [stream_axis,close] = plt.subplots(1, 2)

    ballistic_data = single_ballistic_data

    radial_graph = np.array(ballistic_data["radial_graph"])
    stream_data = np.array(ballistic_data["stream_height"])
    stream_axis.plot(radial_graph, stream_data,color="orange",label="Debris Stream Width")
    stream_axis.set_xlabel("radial coordinate")
    stream_axis.set_ylabel("stream width")


    close.plot(radial_graph, stream_data,color="orange",label="Debris Stream Width")
    close.set_xlabel("radial coordinate")
    close.set_ylabel("stream width")
    close.set_title("Close up view")
    close.set_xlim([15,25])
    close.set_ylim([0,0.14])



    fig.suptitle("Debris Stream Width Evolution")
    plt.tight_layout()
    plt.legend()
    plt.savefig("graphs/stream_height")

def compare_components(a, start_search_divisions):
    fig, axs = plt.subplots(2, 4, figsize=(9, 5))

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

def compare_single_geodesic_components(a, start_search_divisions):
    fig, axs = plt.subplots(2, 4, figsize=(9, 5))


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



def plot_components():
    fig, axs = plt.subplots(1, 4, figsize=(9, 5))

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


def plot_geodesic_components(mino):
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

def plot_3d_star_chunks():
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
            colors.append((index/len(star_chunks),0,1))
            color_list.append(star_chunk["binding_energy"])
 # Red, Green, Blue
    n_bins = 100  # Number of color gradations
    custom_cmap = LinearSegmentedColormap.from_list("custom", colors, N=n_bins)
    scatter = ax.scatter(x_ends,y_ends,z_ends,c=color_list,cmap=custom_cmap)
    cbar = plt.colorbar(scatter,fraction=0.026, pad=0.08)

    cbar.set_label("Specific Binding Energy")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title("Angular Debris Spread due to Varying Binding Energies")

    plt.tight_layout()
    plt.savefig("graphs/three_d_trajectories")


def plot_3d_geodesic():



    ax = plt.axes(projection='3d')

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

    ax.plot(x,y,z,color=colors["blue"],label="Center-of-mass geodesic")
    ax.set_title("Elucidated 3-dimensional Orbital Trajectory",size="large")
    plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')

    #plt.legend(loc="upper right")




    plt.savefig("graphs/single_three_d_trajectory")
    fig, ax2 = plt.subplots(1, 1, figsize=(3, 3))


    ax2.plot(y,z,color=colors["blue"],label="Projected center-of-mass geodesic")
    ax2.scatter(0,0,color="black",label="SMBH location")

    ax2.set_xlim([-0.2,0.2])
    ax2.set_ylim([-12,12])
    ax2.set_ylabel("z",size="xx-large")
    ax2.set_xlabel("x",size="xx-large")
   # plt.legend()
    plt.tight_layout()
   # ax2.set_title("Close-up Side View of the Stream Precession")


    plt.savefig("graphs/single_three_d_trajectory_projection")


def plot_stellar_profile():
    fraction_masses = []
    binding_energies = []

    for star_chunk in star_chunks:
        fraction_of_star = star_chunk["fraction_of_star"]
        binding_energy = star_chunk["binding_energy"]
        fraction_masses.append(fraction_of_star * len(star_chunks))
        binding_energies.append(float(binding_energy))

    stellar_profile, (stellar_profile_axis) = plt.subplots(1, 1, figsize=(6.4, 4), layout="constrained")
    stellar_profile_axis.xaxis.set_inverted(True)
    stellar_profile_axis.plot(binding_energies, fraction_masses)
    stellar_profile.suptitle("Fraction of Stellar Mass at Each Binding Energy")
    stellar_profile_axis.set_xlabel('Binding Energy of Slices')
    plt.xticks(np.linspace(binding_energies[0],binding_energies[-1],5))
    stellar_profile_axis.set_ylabel('Fraction of Stellar Mass')
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

def plot_angular_momentum_over_time():
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
    axs[0].set_xlabel("time")
    axs[0].set_ylabel("Specific Angular Momentum")
    axs[0].plot(times,xycomps,color=colors["green"], label="xy-component of total angular momentum")
    axs[1].plot(times,angles,color=colors["orange"],label="Cosine of approximate total debris inclination")
    axs[1].set_xlabel("time")
    axs[1].set_ylabel("$L_z/L$")

    fig.legend(loc="right")
    fig.suptitle("Time Evolution of Total Angular Momentum before Initial Self Intersection")
    plt.tight_layout()

plot_3d_geodesic()


plt.show()
