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

stations = np.loadtxt(f'{KAYAK_NAME}Stations.csv', delimiter=',')
chinesY = np.loadtxt(f'{KAYAK_NAME}HB.csv', delimiter=',')
chinesZ = np.loadtxt(f'{KAYAK_NAME}HAB.csv', delimiter=',')
keelZ = np.loadtxt(f'{KAYAK_NAME}Keel.csv', delimiter=',')
deckridgeZ = np.loadtxt(f'{KAYAK_NAME}Deckridge.csv', delimiter=',')
gunwaleY = np.loadtxt(f'{KAYAK_NAME}GunwaleHB.csv', delimiter=',')
gunwaleZ = np.loadtxt(f'{KAYAK_NAME}GunwaleHAB.csv', delimiter=',')

if chinesY.ndim == 1:
    chinesY.shape = (1, len(chinesY))
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

def addPointsToDocument(points, feature_name):
    ptList = []
    for i, pt in enumerate(points):
        vec = FreeCAD.Vector(pt[1] * 10, pt[0] * 10, pt[2] * 10)
        obj = doc.addObject("Part::Feature", f"{feature_name}_{i+1}")
        obj.Shape = Part.Vertex(vec)
        obj.Placement = FreeCAD.Placement(vec, FreeCAD.Rotation(0,0,0,1))
        obj.Visibility = True
        ptList.append(obj)
    return ptList

def interpolateSplineFromPoints(points, feature_name):
    spline = Part.BSplineCurve()
    spline.interpolate([pt.Shape.Vertexes[0].Point for pt in points])
    spline_obj = doc.addObject("Part::Feature", feature_name)
    spline_obj.Shape = spline.toShape()
    spline_obj.Visibility = True
    return spline_obj

def twoDOffsetSpline(spline_obj, offset_distance, feature_name):
    offset_obj = FreeCAD.ActiveDocument.addObject("Part::Offset2D", feature_name)
    offset_obj.Source = spline_obj
    offset_obj.Value = offset_distance
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

doc = FreeCAD.newDocument(f'{KAYAK_NAME}')

YZplane = skso.Plane([0,0,0], normal=[1,0,0])
XZplane = skso.Plane([0,0,0], normal=[0,1,0])
XYplane = skso.Plane([0,0,0], normal=[0,0,1])

for idx in range(len(chinesY)):
    chineY = chinesY[idx]
    chineZ = chinesZ[idx]

    points = skso.Points([[stations[i], chineY[i], chineZ[i]] for i in range(len(stations))])

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

    ##### Use the physical model to model the chine as well
    #phys_model = wood_shape.BentWoodShape(twodpoints, tension=1.0, resolution=100)

    ##### Project the physical model into 3d space
    #phys_model_3d = skso.Points([local_to_global(pt, *local_coords) for pt in phys_model.get_curve_points()])

    ##### Find the chine arc endpoints
    #center = local_to_global(local_center, *local_coords)
    #ix1, ix2 = circle.intersect_line(skso.Line([0,0], [1,0]))

    ##### add 3D points to the FreeCAD document
    pts = addPointsToDocument(threedpoints, f'Chine{idx+1}')
    spline = interpolateSplineFromPoints(pts, f'Chine{idx+1}_Spline')
    offset = twoDOffsetSpline(spline, offset_distance=2*25.4, feature_name=f'Chine{idx+1}_Offset')
    extrude = extrudeProfile(offset, length=0.5*25.4, feature_name=f'Chine{idx+1}_Extrude', symmetric=False)


gunwalePoints = skso.Points([[stations[i], gunwaleY[i], gunwaleZ[i]] for i in range(len(stations))])
gunwale_plane = skso.Plane.best_fit(gunwalePoints)
gunwalePointsProjected = skso.Points([gunwale_plane.project_point(pt) for pt in gunwalePoints])
pts = addPointsToDocument(gunwalePointsProjected, 'Gunwale')
spline = interpolateSplineFromPoints(pts, 'Gunwale_Spline')
offset = twoDOffsetSpline(spline, offset_distance=(2*25.4), feature_name='Gunwale_Offset')
extrude = extrudeProfile(offset, length=(0.5*25.4), feature_name='Gunwale_Extrude', symmetric=False, reversed=True)

print(f"Gunwale\n{'-' * 20}")
printPointDistances(gunwalePoints, gunwalePointsProjected)

keelPoints = skso.Points([[stations[i], 0, keelZ[i]] for i in range(len(stations))])
pts = addPointsToDocument(keelPoints, 'Keel')
spline = interpolateSplineFromPoints(pts, 'Keel_Spline')
offset = twoDOffsetSpline(spline, offset_distance=2*25.4, feature_name='Keel_Offset')
extrude = extrudeProfile(offset, length=0.5*25.4, feature_name='Keel_Extrude', symmetric=True)

deckridgePoints = skso.Points([[stations[i], 0, deckridgeZ[i]] for i in range(len(stations))])
addPointsToDocument(deckridgePoints, 'Deckridge')

doc.recompute()
doc.saveAs(f'{KAYAK_NAME}.FCStd')