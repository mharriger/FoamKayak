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
import math
from TableOfOffsets import TableOfOffsets
import skspatial.objects as skso
from skspatial.plotting import plot_3d

stations = np.loadtxt('SeaBeeStations.csv', delimiter=',')
chinesY = np.loadtxt('SeaBeeHB.csv', delimiter=',')
chinesZ = np.loadtxt('SeaBeeHAB.csv', delimiter=',')

offsets = TableOfOffsets(
    stations=list(stations),
    keel=list((None,) * len(stations)),
    chines=[[[chinesY[i][j], chinesZ[i][j]] for j in range(len(chinesY[0]))] for i in range(len(chinesY))],
    gunwale=list((None,) * len(stations)),
    deckridge=list((None,) * len(stations))
    )

fig2d = plt.figure()
fig3d = plt.figure()

axs2d = fig2d.subplots(nrows=len(offsets.chines), ncols=1, sharex=True,
                         width_ratios=[1])

axs3d = fig3d.subplots(nrows=1, ncols=1, subplot_kw={'projection': '3d'})
fig3d.tight_layout()

def global_to_local(point, origin, u, v):
    return np.dot((point - origin), u), np.dot((point - origin), v)

def local_to_global(point, origin, u, v):
    return origin + point[0] * u + point[1] * v

def make_coordinate_system(point, normal):
    chine_plane.normal.unit = normal / np.linalg.norm(normal)
    d = -1 * (normal[0] * point[0] + normal[1] * point[1] + normal[2] * point[2])
    zo = d/(-1 * normal[2])
    print (zo)
    po = np.array((0, 0, zo))

    # Another point on the X axis, x=100, y = 0
    z100 = ((normal[0] * 100) + d)/(-1 * normal[2])
    #z100 = (100 * best_eq[0] + best_eq[2]) / best_eq[3]

    #Local x axis
    xl = np.array([100, 0, z100]) - po
    xl = xl / np.linalg.norm(xl) 

    #Local y axis
    yl = np.cross(xl, chine_plane.normal.unit)
    yl = yl / np.linalg.norm(yl)

    return po, xl, yl

for idx in range(len(chinesY)):
    chineY = chinesY[idx]
    chineZ = chinesZ[idx]

    points = skso.Points([[stations[i], chineY[i], chineZ[i]] for i in range(len(stations))])

    ##### Begin plane fitting

    chine_plane = skso.Plane.best_fit(points)

    ##### End plane fitting

    ##### Create a local coordinate system for the plane
    po, xl, yl = make_coordinate_system(chine_plane.point, chine_plane.normal)

    ##### Convert all of the points on this chine to the LCS
    twodstations = [chine_plane.intersect_plane(skso.Plane(point=[x,0,0], normal=[1,0,0])) for x in stations]
    threedpoints = chine_plane.project_points(points)
    twodpoints = skso.Points([global_to_local(point, po, xl, yl) for point in threedpoints])

    circle = skso.Circle.best_fit(twodpoints)
    local_center = circle.point
    radius = circle.radius

    points_on_circle = skso.Points([circle.project_point(pt) for pt in twodpoints])
    points_on_circle_3d = skso.Points([local_to_global(pt, po, xl, yl) for pt in points_on_circle])

    center = local_to_global(local_center, po, xl, yl)
    ix1, ix2 = circle.intersect_line(skso.Line([0,0], [1,0]))

    ##### 2D plot of chine
    ax = axs2d[idx]
    ax.set_box_aspect(1/10)
    ax.set_ylim((0, 40))
    ax.set_xlim(-10, 400)
    ax.set_xlabel('$x$ (cm)',fontsize=16)
    ax.set_ylabel('\n$y$ (cm)',fontsize=16)
    ax.grid(True, which='both')
    ax.set_aspect('equal')

    twodpoints.plot_2d(ax, c='r')
    points_on_circle.plot_2d(ax, c='g')
    for station in twodstations:
        station.plot_2d(ax, c='g', t_1=-20, t_2=20)

    #draw arc
    arc_angles = np.linspace(np.arccos((ix1 - local_center[0])/radius), np.arccos((ix2 - local_center[0])/radius), 2000)
    arc_xs = (radius * np.cos(arc_angles)) + local_center[0]
    arc_ys = (radius * np.sin(arc_angles)) - abs(local_center[1])
    ax.plot(arc_xs, arc_ys, color = 'c', lw = 3)

    #3D plot of kayak
    points.plot_3d(axs3d, c='b')
    threedpoints.plot_3d(axs3d, c='r')
    points_on_circle_3d.plot_3d(axs3d, c='g')

    axs3d.set_zlim3d(bottom=0, top=20)
    axs3d.set_ylim3d(bottom=0, top=20)
    #axs3d.set_xlim3d(left=50, right=70)

    print("-" * 10, "Chine", idx, "-" * 10)
    print(chine_plane.normal)
    print(points)
    print(threedpoints)
    print (points - points_on_circle_3d)

axs3d.set_aspect('equal')
plt.show()


#   3D plot 
#fig = plt.figure()
#ax = fig.add_subplot()
#ax.plot_surface(xx, yy, z)
#ax.scatter(correctX[1:-1], correctY, correctZ, zdir='z', s=20, c='b',rasterized=True)
#ax.set_xlim3d(-100,1000)
#ax.set_ylim3d(-100, 1000)
#ax.set_zlim3d(-100, 1000)
#ax.set_xlabel('$x$ (cm)',fontsize=16)
#ax.set_ylabel('\n$y$ (cm)',fontsize=16)
#zlabel = ax.set_zlabel('\n$z$ (mm)',fontsize=16)

#ax2 = fig.add_subplot()
#arr = np.array((p0l, p1l, p2l)).T
#print(arr)
#ax2.plot(*arr, 'ro')
#ax2.set_aspect('equal', 'box')
plt.show()
# plt.savefig('steelBallFitted.pdf', format='pdf', dpi=300, bbox_extra_artists=[zlabel], bbox_inches='tight')