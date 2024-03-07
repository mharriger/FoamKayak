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
from mpl_toolkits.mplot3d import Axes3D
import math
from TableOfOffsets import TableOfOffsets
import pyransac3d as pyrsc
import Geometry3D

stations = np.loadtxt('SeaBeeStations.csv', delimiter=',')
chinesY = np.loadtxt('SeaBeeHB.csv', delimiter=',')
chinesZ = np.loadtxt('SeaBeeHAB.csv', delimiter=',')

#chine3d = np.array((stations[1:-1], chinesY, chinesZ))
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

axsBottom = fig3d.subplots(nrows=1, ncols=1, subplot_kw={'projection': '3d'})
fig3d.tight_layout()

for frame in offsets.frames:
    for centroid in frame:
        print(centroid)
        axsBottom.plot(*centroid, 'bo')


# https://stackoverflow.com/questions/28910718/give-3-points-and-a-plot-circle
def define_circle(p1, p2, p3):
    """
    Returns the center and radius of the circle passing the given 3 points.
    In case the 3 points form a line, returns (None, infinity).
    """
    temp = p2[0] * p2[0] + p2[1] * p2[1]
    bc = (p1[0] * p1[0] + p1[1] * p1[1] - temp) / 2
    cd = (temp - p3[0] * p3[0] - p3[1] * p3[1]) / 2
    det = (p1[0] - p2[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p2[1])
    
    if abs(det) < 1.0e-6:
        return (None, np.inf)
    
    # Center of circle
    cx = (bc*(p2[1] - p3[1]) - cd*(p1[1] - p2[1])) / det
    cy = ((p1[0] - p2[0]) * cd - (p2[0] - p3[0]) * bc) / det
    
    radius = np.sqrt((cx - p1[0])**2 + (cy - p1[1])**2)
    return ((cx, cy), radius)

def plane3points(p0, p1, p2):
    x0, y0, z0 = p0
    x1, y1, z1 = p1
    x2, y2, z2 = p2

    ux, uy, uz = [x1-x0, y1-y0, z1-z0]
    vx, vy, vz = [x2-x0, y2-y0, z2-z0]

    u_cross_v = [uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx]

    point = np.array((x0, y0, z0))

    normal = np.array(u_cross_v)
    return point, normal

def global_to_local(point, origin, u, v):
    return np.dot((point - origin), u), np.dot((point - origin), v)

def local_to_global(point, origin, u, v):
    return origin + point[0] * u + point[1] * v

def intersect_planes(point, normal1, normal2):
    v = normal1 * normal2
    return (point, point + v)

def intersect_x_axis_circle(c, r):
    return c[0] - math.sqrt(r**2 - c[1]**2), c[0] + math.sqrt(r**2 - c[1]**2)

def make_coordinate_system(point, normal):
    unit_normal = normal / np.linalg.norm(normal)
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
    yl = np.cross(xl, unit_normal)
    yl = yl / np.linalg.norm(yl)

    return po, xl, yl

def intersect_plane_circle(c, r, c_normal, p, p_normal):
    # Make a local coordinate system (LCS)
    po, xl, yl = make_coordinate_system(c, c_normal)



for idx in range(len(chinesY)):
    chineY = chinesY[idx]
    chineZ = chinesZ[idx]
    print(chineY, chineZ)
    
    points = np.array([[stations[i], chineY[i], chineZ[i]] for i in range(len(stations))])

    ##### Begin plane fitting. TODO: A more accurate best fit
    # Pick 3 points
    p0 = np.array((stations[0], chineY[0], chineZ[0]))
    p1 = np.array((stations[2], chineY[2], chineZ[2]))
    p2 = np.array((stations[4], chineY[4], chineZ[4]))

    circle1 = pyrsc.Circle()
    center, axis, radius, inliers = circle1.fit(points, 1)
    normal = axis

    print('RANSAC results: ', axis, radius, center, inliers)

#    # Define the plane containing those points
    point, normal_old = plane3points(p0, p1, p2)
#    unit_normal = normal / np.linalg.norm(normal)
    unit_normal = axis / np.linalg.norm(axis)
    unit_normal_old = normal_old / np.linalg.norm(normal_old)

    chine_plane = Geometry3D.Plane(Geometry3D.Point(center[0], center[1], center[2]), Geometry3D.Vector(unit_normal[0], unit_normal[1], unit_normal[2]))
    print ("chine plane parameters:", chine_plane.parametric())
#    ##### End plane fitting
#
    print(unit_normal, unit_normal_old)
#
    ##### Create a local coordinate system for the plane
    po, xl, yl = make_coordinate_system(center, axis)

    local_center = global_to_local(center, po, xl, yl)
    ix1, ix2 = intersect_x_axis_circle(local_center, radius)

    ##### Convert all of the points on this chine to the LCS
    #TODO: Find nearest point on plane to actual point. Not sure what this actually gives although it looks close
    twodpoints = []
    for i in range(len(chineY)):
        centroid = np.array((stations[i],chineY[i],chineZ[i]))
        twodpoints.append(global_to_local(centroid, po, xl, yl))

    ### Find circle through the three translated points
    c, r = define_circle(twodpoints[0], twodpoints[math.ceil(len(twodpoints) / 2.0)], twodpoints[-1])
    print("c,r:",c,r)

    # Convert c to global coords for display
    C = local_to_global(c, po, xl, yl)
    print("global center: ", C)

    # Determine station locations as lines on the chine plane
    for x in stations:
        station_plane = Geometry3D.Plane(Geometry3D.Point(x, 0, 0), Geometry3D.Vector(1, 0, 0))
        station_line = station_plane.intersection(chine_plane)
        so = station_line.intersection(Geometry3D.Plane.xz_plane())
        print("Station Origin ", so)
        print("Station Origin Local ", global_to_local(np.array((so.x, so.y, so.z)), po, xl, yl))
        #TODO: How to intersect line with circle (do it in LCS?). Geometry3D circle is actually a bunch of segments.
        # Also, how to plot chine line on chine plane diagram
    
    nyz = np.array((1, 0, 0))
    station_vec = nyz * unit_normal
    print("station_vec:", station_vec)

#    twodstations = []
#    for x in stations[1:-1]:
#        zl = (normal[0] * x + d)/(-normal[2])
#        twodstations.append(np.array((
#            global_to_local(po, np.array((-x, 0, zl)), xl, yl),
#            global_to_local(po, np.array((-x, 0, zl)) - station_vec, xl, yl) * np.array((1, -100))
#        )))
#    print(twodstations)
    # Adjust other chine points to intersection of station plane with chine arc

    # Calculate starting point and angle of arc
    # Starts and ends at points circle intersects local X axis
    print("ix1, ix2:",ix1,ix2)

    c0 = np.array((ix1, 0))


    ##### 2D plot of chine
    ax = axs2d[idx]
    ax.set_box_aspect(1/10)
    ax.set_ylim((0, 40))
#    ax.plot(*zip(*twodpoints), 'ro')
    #ax.plot(*zip(p0a, p1a, p2a), 'ro')
#    for station in twodstations:
#        ax.plot(*zip(*station), 'g')
    ax.plot([ix1, ix2], [0, 0], 'ro')
    ax.set_xlim(-10, 400)
    #ax.set_ylim(-5, 40)
    ax.set_xlabel('$x$ (cm)',fontsize=16)
    ax.set_ylabel('\n$y$ (cm)',fontsize=16)
    ax.grid(True, which='both')
    ax.set_aspect('equal')

    #draw arc
    arc_angles = np.linspace(np.arccos((ix1 - local_center[0])/radius), np.arccos((ix2 - local_center[0])/radius), 2000)
    arc_xs = (radius * np.cos(arc_angles)) + local_center[0]
    arc_ys = (radius * np.sin(arc_angles)) - abs(local_center[1])
    ax.plot(arc_xs, arc_ys, color = 'c', lw = 3)

#    #3D
#    for centroid in twodpoints:
#        print(local_to_global(centroid, po, xl, yl))
#        axsRight.plot(*local_to_global(centroid, po, xl, yl), 'ro')
#    axsRight.plot(*local_to_global((ix1, 0), po, xl, yl), 'ro')
#    axsRight.plot(*local_to_global((ix2, 0), po, xl, yl), 'ro')
#
#
    for station_x in offsets._stations:
        station_plane = plane3points(np.array((station_x, 0, 0)), np.array((station_x, 100, 0)), np.array((station_x, 0, 100)))



axsBottom.set_aspect('equal')
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