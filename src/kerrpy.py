import kerrgeopy as kg
from math import cos, pi
import numpy as np
import matplotlib.pyplot as plt
#orbit = kg.StableOrbit.from_constants(0.9, 0.9, 4.3, 2.5)
orbit = kg.StableOrbit.from_constants(0.9, 0.95, 1.6, 8)
print(orbit.constants_of_motion())
t, r, theta, phi = orbit.trajectory()


fig, ax = orbit.plot(0,10)


time = np.linspace(0,10,200)

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





plt.show()