# Goal is:
#   Input: Table of offsets in Stations/HB/HAB
#   Output:
#       Updated offsets
#       3D model of all components
#       Scaled PDF diagrams of components with measurements

import sys
sys.path.append("/usr/lib/freecad-python3/lib")
sys.path.append("/usr/share/freecad/Mod")

import FreeCAD
import FreeCADGui as Gui
import PartDesign
import Part


import numpy as np
#from TableOfOffsets import TableOfOffsets
import skspatial.objects as skso

KAYAK_NAME = 'SeaBee'

YZplane = skso.Plane([0,0,0], normal=[1,0,0])
XZplane = skso.Plane([0,0,0], normal=[0,1,0])
XYplane = skso.Plane([0,0,0], normal=[0,0,1])

stations = np.loadtxt(f'{KAYAK_NAME}Stations.csv', delimiter=',')
chinesX = np.loadtxt(f'{KAYAK_NAME}HB.csv', delimiter=',')
chinesZ = np.loadtxt(f'{KAYAK_NAME}HAB.csv', delimiter=',')
keelZ = np.loadtxt(f'{KAYAK_NAME}Keel.csv', delimiter=',')
deckridgeZ = np.loadtxt(f'{KAYAK_NAME}Deckridge.csv', delimiter=',')
gunwaleX = np.loadtxt(f'{KAYAK_NAME}GunwaleHB.csv', delimiter=',')
gunwaleZ = np.loadtxt(f'{KAYAK_NAME}GunwaleHAB.csv', delimiter=',')

if chinesX.ndim == 1:
    chinesX.shape = (1, len(chinesX))
if chinesZ.ndim == 1:
    chinesZ.shape = (1, len(chinesZ))
#
#offsets = TableOfOffsets(
#    stations=list(stations),
#    keel=list((None,) * len(stations)),
#    chines=[[[chinesY[i][j], chinesZ[i][j]] for j in range(len(chinesY[0]))] for i in range(len(chinesY))],
#    gunwale=list((None,) * len(stations)),
#    deckridge=list((None,) * len(stations))
#    )
#
def global_to_local(point: skso.Point, origin: skso.Point, u: skso.Vector, v: skso.Vector) -> skso.Point:
    # Ensure numeric arrays and return 2D coordinates in the plane basis (u,v)
    p = np.asarray(point, dtype=float) - np.asarray(origin, dtype=float)
    u_arr = np.asarray(u, dtype=float)
    v_arr = np.asarray(v, dtype=float)
    return skso.Point([float(np.dot(p, u_arr)), float(np.dot(p, v_arr))])

def local_to_global(point: skso.Point, origin: skso.Point, u: skso.Vector, v: skso.Vector) -> skso.Point:
    # Accept a 2D local point and return a 3D point on the plane
    local = np.asarray(point, dtype=float)
    origin_arr = np.asarray(origin, dtype=float)
    u_arr = np.asarray(u, dtype=float)
    v_arr = np.asarray(v, dtype=float)
    global_arr = origin_arr + float(local[0]) * u_arr + float(local[1]) * v_arr
    return skso.Point(global_arr)

def make_coordinate_system(plane: skso.Plane) -> tuple[skso.Point, skso.Vector, skso.Vector]:
    # Build a stable orthonormal basis on the plane.
    # Origin: closest point on the plane to world origin (projection of origin onto plane)
    origin = plane.project_point(skso.Point([0.0, 0.0, 0.0]))

    # Normal vector (unit)
    n = np.asarray(plane.normal, dtype=float)
    n = n / np.linalg.norm(n)

    # Choose a reference axis that's not parallel to the normal
    ref = np.array([0.0, 0.0, 1.0], dtype=float)
    if abs(np.dot(n, ref)) > 0.999:
        ref = np.array([0.0, 1.0, 0.0], dtype=float)

    # u is one axis in the plane, v is the other; ensure right-handed basis: u x v = n
    u = np.cross(ref, n)
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    v /= np.linalg.norm(v)

    return skso.Point(origin), skso.Vector(u), skso.Vector(v)

def addPointsToDocument(points, feature_name):
    ptList = []
    for i, pt in enumerate(points):
        vec = FreeCAD.Vector(pt[0] , pt[1] , pt[2] )
        obj = doc.addObject("Part::Vertex", f"{feature_name}_{i+1}")
        obj.Shape = Part.Vertex(vec)
        obj.Placement = FreeCAD.Placement(vec, FreeCAD.Rotation(0,0,0,1))
        obj.Visibility = True
        ptList.append(obj)
    return ptList

def interpolateSplineFromPoints(points, feature_name):
    spline = Part.BSplineCurve()
    spline.interpolate([pt.Shape.Vertexes[0].Point for pt in points])
    spline_obj = doc.addObject("Part::Part2DObject", feature_name)
    spline_obj.Shape = spline.toShape()
    return spline_obj

def twoDOffsetSpline(spline_obj, offset_distance, feature_name):
    offset_obj = FreeCAD.ActiveDocument.addObject("Part::Offset2D", feature_name)
    offset_obj.Source = spline_obj
    offset_obj.Value = -1.0 * offset_distance
    offset_obj.Mode = "Skin"
    offset_obj.Fill = True
    return offset_obj

def extrudeProfile(profile_obj, length, feature_name, symmetric=False, reversed=False):
    extrude_obj = FreeCAD.ActiveDocument.addObject("Part::Extrusion", feature_name)
    extrude_obj.Base = profile_obj
    extrude_obj.DirMode = "Normal"
    extrude_obj.LengthFwd = length
    extrude_obj.Symmetric = symmetric
    extrude_obj.Reversed = reversed
    return extrude_obj

def printPointDistances(points1, points2):
    if len(points1) != len(points2):
        raise ValueError("Point lists must be of the same length")
    for i in range(len(points1)):
        print(skso.Point(points1[i]).distance_point(skso.Point(points2[i])))

def find_chine_endpoints(twodpoints, chine_plane, local_coords):
    # Fit a circle in 2D, create corresponding 3D sphere and intersect with plane intersection line
    circle = skso.Circle.best_fit(twodpoints)
    distances = [circle.distance_point(pt) for pt in twodpoints]
    line = chine_plane.intersect_plane(YZplane)
    circle_3d_center = local_to_global(circle.point, *local_coords)
    sphere = skso.Sphere(circle_3d_center, circle.radius)
    endpoints = sphere.intersect_line(line)
    if len(endpoints) < 2:
        raise RuntimeError("Could not find two intersection points between sphere and line")
    return distances, endpoints

doc = FreeCAD.newDocument(f'{KAYAK_NAME}')

for idx in range(len(chinesX)):
    chineX = chinesX[idx]
    chineZ = chinesZ[idx]

    points = skso.Points([[chineX[i], stations[i], chineZ[i]] for i in range(len(stations))])

    ##### Find the best-fit plane
    chine_plane = skso.Plane.best_fit(points)

    ##### Create a local coordinate system for the plane
    local_coords = make_coordinate_system(chine_plane)

    ##### Project all of the points on this chine to the plane
    threedpoints = [chine_plane.project_point(pt) for pt in points]
    twodpoints = skso.Points([global_to_local(skso.Point(point), *local_coords) for point in threedpoints])

    print(f"Chine {idx+1}\n{'-' * 20}")
    printPointDistances(points, threedpoints)
    print("-" * 20)

    distances, chine_endpoints = find_chine_endpoints(twodpoints, chine_plane, local_coords)
    print(f"Circle fit distances: {distances}")
    threedpoints.insert(0, chine_endpoints[0])
    threedpoints.append(chine_endpoints[1])

    ##### Use the physical model to model the chine as well
    #phys_model = wood_shape.BentWoodShape(twodpoints, tension=1.0, resolution=100)

    ##### Project the physical model into 3d space
    #phys_model_3d = skso.Points([local_to_global(pt, *local_coords) for pt in phys_model.get_curve_points()])

    ##### Find the chine arc endpoints
    #center = local_to_global(local_center, *local_coords)
    #ix1, ix2 = circle.intersect_line(skso.Line([0,0], [1,0]))

    ##### add 3D points to the FreeCAD document
    pts = addPointsToDocument(threedpoints, f'Chine{idx+1}')
    doc.recompute()
    spline = interpolateSplineFromPoints(pts, f'Chine{idx+1}_Spline')
    offset = twoDOffsetSpline(spline, offset_distance=2*25.4, feature_name=f'Chine{idx+1}_Offset')
    extrude = extrudeProfile(offset, length=0.5*25.4, feature_name=f'Chine{idx+1}_Extrude', symmetric=False)


gunwalePoints = skso.Points([[gunwaleX[i], stations[i], gunwaleZ[i]] for i in range(len(stations))])
gunwale_plane = skso.Plane.best_fit(gunwalePoints)
gunwalePointsProjected = [gunwale_plane.project_point(pt) for pt in gunwalePoints]
print(f"Gunwale\n{'-' * 20}")
printPointDistances(gunwalePoints, gunwalePointsProjected)
local_coords = make_coordinate_system(gunwale_plane)
gunwalePoints2d = [global_to_local(skso.Point(point), *local_coords) for point in gunwalePointsProjected]
distances, endpoints = find_chine_endpoints(gunwalePoints2d, gunwale_plane, local_coords)
gunwalePointsProjected.insert(0, endpoints[0])
gunwalePointsProjected.append(endpoints[1])
pts = addPointsToDocument(gunwalePointsProjected, 'Gunwale')
doc.recompute()
spline = interpolateSplineFromPoints(pts, 'Gunwale_Spline')
offset = twoDOffsetSpline(spline, offset_distance=(2*25.4), feature_name='Gunwale_Offset')
extrude = extrudeProfile(offset, length=(0.5*25.4), feature_name='Gunwale_Extrude', symmetric=False, reversed=True)

keelPoints = skso.Points([[0, stations[i], keelZ[i]] for i in range(len(stations))])
pts = addPointsToDocument(keelPoints, 'Keel')
doc.recompute()
spline = interpolateSplineFromPoints(pts, 'Keel_Spline')
offset = twoDOffsetSpline(spline, offset_distance=2*25.4, feature_name='Keel_Offset')
extrude = extrudeProfile(offset, length=0.5*25.4, feature_name='Keel_Extrude', symmetric=True)

deckridgePoints = skso.Points([[0, stations[i], deckridgeZ[i]] for i in range(len(stations))])
addPointsToDocument(deckridgePoints, 'Deckridge')

#TODO: Make a plane for each station
for station in stations:
    station_plane = FreeCAD.ActiveDocument.addObject("Part::Plane", f"StationPlane_{station}")
    station_plane.Length = 2000  # large enough to cover the model
    station_plane.Width = 1500
    station_plane.Placement = FreeCAD.Placement(FreeCAD.Vector(-1000, station, -500), FreeCAD.Rotation(1,0,0,1))

doc.recompute()
doc.saveAs(f'{KAYAK_NAME}.FCStd')