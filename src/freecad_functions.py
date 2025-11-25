# Function to interact with FreeCAD
import numpy as np

import sys
sys.path.append("/usr/lib/freecad-python3/lib")
sys.path.append("/usr/share/freecad/Mod")
import FreeCAD as App
import FreeCADGui as Gui
import Part
import Sketcher

from Bspline import Bspline

def makeFreeCADDocument(name):
    doc = App.newDocument(name)
    return doc

def add_bspline_sketch(doc, name, point, normal, bspline: Bspline, chine_points_2d, rotation=None):

    conList = []
    sketch = doc.addObject('Sketcher::SketchObject', name)
    # Set the sketch plane
    normal = normal / np.linalg.norm(normal)
    if rotation is not None:
        rotation = App.Rotation(App.Vector(rotation[0], rotation[1], rotation[2]), rotation[3])
    else:
        fc_pl = Part.Plane(App.Vector(*point), App.Vector(*normal))
        
        ###DEBUG###
        o = doc.addObject('Part::Plane', f'{name}_Plane')
        o.Length = 2000  # large enough to cover the model
        o.Width = 1500
        o.Placement = App.Placement(App.Vector(*point), fc_pl.Rotation) 
        ###END_DEBUG###
        rotation = fc_pl.Rotation
    sketch.Placement = App.Placement(App.Vector(*point), rotation)
    # Add the B-spline to the sketch
    points_2d = [App.Vector(abs(pt[1]), abs(pt[0]), 0) for pt in bspline.control_points]
    curve = Part.BSplineCurve()
    curve.increaseDegree(3)
    curve.buildFromPolesMultsKnots(points_2d, bspline.multiplicities, bspline.knots, False, bspline.degree)
    curve_n = sketch.addGeometry(curve)
    for pt_idx, pt in enumerate(chine_points_2d):
        y,x = [abs(pt[0]), abs(pt[1])]
        n = sketch.addGeometry(Part.Point(App.Vector(x,y)),True)
        conList.append(Sketcher.Constraint('InternalAlignment:Sketcher::BSplineKnotPoint',n,1,curve_n,pt_idx))
        conList.append(Sketcher.Constraint('DistanceX',-1,1,n,1,x))
        conList.append(Sketcher.Constraint('DistanceY',-1,1,n,1,y))
    sketch.addConstraint(conList)
    conList = []
    sketch.exposeInternalGeometry(curve_n)
    # First and last control points coincident with first and last b-spline knots
    controlPoints = filter(lambda g: g.TypeId == 'Part::GeomCircle', sketch.Geometry)
    id = next(controlPoints).getExtensionOfType('Sketcher::SketchGeometryExtension').Id - 1
    sketch.addConstraint(Sketcher.Constraint('Coincident', id + 1, 3, 0, 1))  # First control point
    id = list(controlPoints)[-1].getExtensionOfType('Sketcher::SketchGeometryExtension').Id - 1
    sketch.addConstraint(Sketcher.Constraint('Coincident', id - 1, 3, 0, 2))  # Last control point
    return sketch