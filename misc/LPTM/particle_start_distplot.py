import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plot
import numpy as np

PARTICLES_FILENAME = "partstartpos"

with open(PARTICLES_FILENAME) as fh:
    particles_data = np.genfromtxt(fh.readlines()[1:], unpack=True)

t_release, x, y, z = particles_data

fig = plot.figure(figsize=(10,4))

plot.subplot(121)
plot.plot(*np.unique(z, return_counts=True), marker='.')
plot.xlabel("height [m]")
plot.ylabel("number of particles [1]")

plot.subplot(122)
xy_pos = zip(x, y)
xy_pos_unique = set(xy_pos)
plot.plot(*xy_pos_unique, marker='.', linestyle='')
plot.xlabel("x-distance [m]")
plot.ylabel("y-distance [m]")
plot.gca().set_aspect(1)

plot.tight_layout()
plot.savefig(__file__.replace('.py', '.png'))
