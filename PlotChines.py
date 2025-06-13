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

import wood_shape

KAYAK_NAME = 'SeaBee'

stations = np.loadtxt(f'{KAYAK_NAME}Stations.csv', delimiter=',')
chinesY = np.loadtxt(f'{KAYAK_NAME}HB.csv', delimiter=',')
chinesZ = np.loadtxt(f'{KAYAK_NAME}HAB.csv', delimiter=',')

if chinesY.ndim == 1:
    chinesY.shape = (1, len(chinesY))
if chinesZ.ndim == 1:
    chinesZ.shape = (1, len(chinesZ))

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

def global_to_local(point: skso.Point, origin: skso.Point, u: skso.Vector, v: skso.Vector) -> skso.Point:
    return skso.Point((np.dot((point - origin), u), np.dot((point - origin), v)))

def local_to_global(point: skso.Point, origin: skso.Point, u: skso.Vector, v: skso.Vector) -> skso.Point:
    return origin + point[0] * u + point[1] * v

def make_coordinate_system(plane: skso.Plane) -> tuple[skso.Point, skso.Vector, skso.Vector]:
    # Origin = intersection of chine plane and z-axis
    po = plane.intersect_line(skso.Line([0,0,0], [0,0,1]))

    # Local x axis = intersection of chine plane and XZ plane
    xl = plane.intersect_plane(skso.Plane([0,0,0], normal=[0,-1,0]))

    #Local y axis
    yl = skso.Vector(np.cross(xl.direction.unit(), plane.normal.unit()))

    return po, xl.direction.unit(), yl.unit()

for idx in range(len(chinesY)):
    chineY = chinesY[idx]
    chineZ = chinesZ[idx]

    points = skso.Points([[stations[i], chineY[i], chineZ[i]] for i in range(len(stations))])

    ##### Find the best-fit plane
    chine_plane = skso.Plane.best_fit(points)

    ##### Create a local coordinate system for the plane
    local_coords = make_coordinate_system(chine_plane)

    ##### Project the station lines and all of the points on this chine to the plane
    twodstations = [chine_plane.intersect_plane(skso.Plane(point=[x,0,0], normal=[1,0,0])) for x in stations]
    threedpoints = chine_plane.project_points(points)
    twodpoints = skso.Points([global_to_local(skso.Point(point), *local_coords) for point in threedpoints])

    ##### Find the circle that best fits the projected points
    circle = skso.Circle.best_fit(twodpoints)
    local_center = circle.point
    radius = circle.radius

    ##### Use the physical model to model the chine as well
    phys_model = wood_shape.BentWoodShape(twodpoints, tension=1.0, resolution=100)

    ##### Project the 2D points to the circle, then reproject back to 3D space
    points_on_circle = skso.Points([circle.project_point(pt) for pt in twodpoints])
    points_on_circle_3d = skso.Points([local_to_global(pt, *local_coords) for pt in points_on_circle])

    ##### Project the physical model into 3d space
    phys_model_3d = skso.Points([local_to_global(pt, *local_coords) for pt in phys_model.get_curve_points()])

    ##### Find the chine arc endpoints
    center = local_to_global(local_center, *local_coords)
    ix1, ix2 = circle.intersect_line(skso.Line([0,0], [1,0]))

    ##### 2D plot of chine
    if chinesY.shape[0] > 1:
        ax = axs2d[idx]
    else:
        ax = axs2d
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

    # Draw the physical model
    skso.Points(phys_model.get_curve_points()).plot_2d(ax, c='m', lw=2, label='Physical Model')

    #3D plot of kayak
    points.plot_3d(axs3d, c='b')
    threedpoints.plot_3d(axs3d, c='r')
    points_on_circle_3d.plot_3d(axs3d, c='g')
    phys_model_3d.plot_3d(axs3d, c='m', lw=2, label='Physical Model')

    axs3d.set_zlim3d(bottom=0, top=20)
    axs3d.set_ylim3d(bottom=0, top=20)
    #axs3d.set_xlim3d(left=50, right=70)

    print("-" * 10, "Chine", idx, "-" * 10)
    print(chine_plane.normal)
    print(points)
    print(points_on_circle_3d)
    print([skso.Point(x[0]).distance_point(skso.Point(x[1])) for x in zip(points, points_on_circle_3d)])
    print("Phys model:" , phys_model.get_curve_points())

    phys_model.draw(f'{KAYAK_NAME}_Chine_{idx}.png', width=800, height=600, margin=50)

axs3d.set_aspect('equal')
plt.show()
