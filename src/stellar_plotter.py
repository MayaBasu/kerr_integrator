import math

import kerrgeopy as kg
import numpy as np
import matplotlib.pyplot as plt
import json
import math as m

with open('../star.json', 'r') as file:
    star_chunks = json.load(file)["star_chunks"]


def stream_evolution():
    fig, stream_axis = plt.subplots(1,1)
    for star_chunk in star_chunks:
        ballistic_data = star_chunk["geodesic_graph"]

        radial_graph = np.array(ballistic_data["radial_graph"])
        stream_data = np.array(ballistic_data["stream_height"])
        stream_axis.plot(radial_graph,stream_data)


def compare_components(a,start_search_divisions):
    fig, axs = plt.subplots(3, 4)
    for star_chunk in star_chunks:
        ballistic_data = star_chunk["geodesic_graph"]
        num_steps = ballistic_data["num_steps"]
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
        r_phase,theta_phase = find_initial_phases(a,E,LZ,C,radial_graph[0],theta_graph[0],start_search_divisions)
        orbit = kg.StableOrbit.from_constants(a, E, LZ, C,initial_phases=(0,r_phase,theta_phase,0))
        t_kg, r_kg, theta_kg, phi_kg = orbit.trajectory()



        x_axis = np.linspace(0, num_steps*step_size, num_steps)

        axs[0,0].plot(x_axis, radial_graph)
        plt.xlabel("$\lambda$")
        plt.ylabel("$r(\lambda)$")

        axs[0,1].plot(x_axis, theta_graph)
        plt.xlabel("$\lambda$")
        plt.ylabel(r"$\theta(\lambda)$")

        axs[0,2].plot(x_axis, phi_graph)
        plt.xlabel("$\lambda$")
        plt.ylabel(r"$\phi(\lambda)$")


        axs[0,3].plot(x_axis, t_graph)
        plt.xlabel("$\lambda$")
        plt.ylabel(r"$t(\lambda)$")


        axs[1,0].plot(x_axis, r_kg(x_axis))
        plt.xlabel("$\lambda$")
        plt.ylabel("$r(\lambda)$")

        axs[1,1].plot(x_axis, theta_kg(x_axis))
        plt.xlabel("$\lambda$")
        plt.ylabel(r"$\theta(\lambda)$")

        axs[1,2].plot(x_axis, phi_kg(x_axis))
        plt.xlabel("$\lambda$")
        plt.ylabel(r"$\phi(\lambda)$")

        axs[1,3].plot(x_axis, t_kg(x_axis))
        plt.xlabel("$\lambda$")
        plt.ylabel(r"$t(\lambda)$")





def plot_3d():
    three_d_fig = plt.figure(figsize=(6, 6))
    ax = plt.axes(projection='3d')
    for star_chunk in star_chunks:
        ballistic_data = star_chunk["geodesic_graph"]

        phi_graph = np.array(ballistic_data["phi_graph"])
        theta_graph = np.array(ballistic_data["theta_graph"])
        radial_graph = np.array(ballistic_data["radial_graph"])

        x = (radial_graph * np.sin(theta_graph) * np.cos(phi_graph))
        y = (radial_graph * np.sin(theta_graph) * np.sin(phi_graph))
        z = (radial_graph * np.cos(theta_graph))
        ax.plot3D(x, y, z)


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
    stellar_profile.suptitle('Fraction of Star with each Binding Energy')
    stellar_profile_axis.plot(binding_energies, fraction_masses)
    stellar_profile.suptitle("Fraction of Stellar Mass at Each Binding Energy")
    stellar_profile_axis.set_xlabel('Binding Energy of Slices')
    stellar_profile_axis.set_ylabel('Fraction of Stellar Mass')
    stellar_profile_axis.grid(True)
    stellar_profile.canvas.manager.set_window_title('fraction_with_binding_energy')
    stellar_profile.tight_layout()
    plt.savefig("graphs/fraction_with_binding_energy.png")



def find_initial_phases(a,E,LZ,C,r0,theta0,start_search_divisions):

    orbit = kg.StableOrbit.from_constants(a, E, LZ, C)
    a,p,e,x = kg.apex_from_constants(a,E,LZ,C)
    t_kg, r_kg, theta_kg, phi_kg = orbit.trajectory()

    half_r_period = math.pi/kg.r_frequency(a,p,e,x)
    half_theta_period = math.pi/kg.theta_frequency(a,p,e,x)
    print(theta0)

    initial_r_index = np.argmin(np.absolute(np.array(r_kg(np.linspace(0,half_r_period,start_search_divisions)))-r0))
    initial_theta_index = np.argmin(np.absolute(np.array(theta_kg(np.linspace(0,half_theta_period,start_search_divisions)))-theta0))

    r_delay_time = initial_r_index*(half_r_period/start_search_divisions)
    theta_delay_time = initial_theta_index*(half_theta_period/start_search_divisions)

    print("r: Matching the initial condition of {} by taking {} at a delay_time of {}, index is {}".format(r0,r_kg(r_delay_time), r_delay_time,initial_r_index))
    print("theta: Matching the initial condition of {} by taking {} at a delay_time of {}, index is {}".format(theta0,theta_kg(theta_delay_time), theta_delay_time,initial_theta_index))

    r_phase = (r_delay_time/half_r_period)*math.pi
    theta_phase = (theta_delay_time/half_theta_period)*math.pi

    return r_phase,theta_phase





compare_components(0.9,100000)


plt.show()
