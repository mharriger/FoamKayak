# Goal is:
#   Input: Table of offsets in Stations/HB/HAB
#   Output:
#       Updated offsets
#       3D model of all components
#       Scaled PDF diagrams of components with measurements

import numpy as np
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
import matplotlib.pyplot as plt
from TableOfOffsets import TableOfOffsets
import skspatial.objects as skso

stations = np.loadtxt('SeaTourStations.csv', delimiter=',')
chinesY = np.loadtxt('SeaTourHB.csv', delimiter=',')
chinesZ = np.loadtxt('SeaTourHAB.csv', delimiter=',')
keel = np.loadtxt('SeaTourKeel.csv', delimiter=',')

if chinesY.ndim == 1:
    chinesY.shape = (1, len(chinesY))
if chinesZ.ndim == 1:
    chinesZ.shape = (1, len(chinesZ))

offsets = TableOfOffsets(
    stations=list(stations),
    keel=list(keel),
    chines=[[[chinesY[i][j], chinesZ[i][j]] for j in range(len(chinesY[0]))] for i in range(len(chinesY))],
    gunwale=list((None,) * len(stations)),
    deckridge=list((None,) * len(stations))
    )

fig2d = plt.figure()
fig3d = plt.figure()

axs2d = fig2d.subplots(nrows=1, ncols=1)

axs3d = fig3d.subplots(nrows=1, ncols=1, subplot_kw={'projection': '3d'})
fig3d.tight_layout()

points = skso.Points(offsets.keel)


##### Find the circle that best fits the keel points
circle = skso.Circle.best_fit(points)

##### 2D plot of keel
axs2d.set_box_aspect(1/10)
axs2d.set_ylim((0, 40))
axs2d.set_xlim(-10, 400)
axs2d.set_xlabel('$x$ (cm)',fontsize=16)
axs2d.set_ylabel('\n$y$ (cm)',fontsize=16)
axs2d.grid(True, which='both')
axs2d.set_aspect('equal')
points.plot_2d(axs2d, c='r')

#draw circle
circle.plot_2d(axs2d)

print([circle.distance_point(skso.Point(p)) for p in points])

plt.show()
