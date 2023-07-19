# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 18:38:32 2023

@author: Ytoga
"""

import numpy as np
from mayavi import mlab

fig = mlab.figure(
    fgcolor = (0, 0, 0),
    bgcolor = (1, 1, 1),
    size = (1000, 1000)
)
fig.scene.parallel_projection = True

u = np.linspace(0, 2*np.pi, 1000)
v = np.linspace(-8, 8, 1000)
x = 2.5 * np.exp(-v**2 / 10) * np.sin(20*u)
y = 2.5 * np.exp(-v**2 / 10) * np.cos(20*u)
z = v

mlab.plot3d(x, y, z, tube_radius=0.09, color=(1, 0, 0))

u, v = np.mgrid[0:np.pi*2:1000j, -8:8:1000j]
x = 2.48 * np.exp(-v**2 / 10) * np.sin(20*u)
y = 2.48 * np.exp(-v**2 / 10) * np.cos(20*u)
z = v
mlab.mesh(x, y, z,
      # color = np.exp(-np.linspace(-8, 8, 1000)**2/10),
      colormap = 'Blues',
      resolution = 20,
      representation = 'surface',
      line_width = 0,
      # tube_radius = 0,
      # tube_sides = 10,
      opacity = 0.9
)
# mlab.savefig('mlaser.png')
mlab.show()
