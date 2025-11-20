import numpy as np
from scipy.interpolate import CubicSpline
import cairo

class BentWoodShape:
    def __init__(self, points, tension=1.0, resolution=100):
        self.points = np.array(points)
        self.tension = tension
        self.resolution = resolution
        self.spline = None
        self.generate_shape()
    
    def generate_shape(self):
        if len(self.points) < 2:
            raise ValueError("At least two points are required")
        
        # Calculate cumulative distance along the curve for parameterization
        t = np.zeros(len(self.points))
        for i in range(1, len(self.points)):
            t[i] = t[i-1] + np.linalg.norm(self.points[i] - self.points[i-1])
        
        # Normalize parameter t
        if t[-1] != 0:
            t = t / t[-1]
        
        # Generate spline with tension
        self.spline = CubicSpline(t, self.points, bc_type='clamped')
        
    def get_curve_points(self):
        t = np.linspace(0, 1, self.resolution)
        return self.spline(t)


def generate_bent_wood_shape(points, tension=1.0, resolution=100):
    """Create a bent wood shape from a list of points"""
    return BentWoodShape(points, tension, resolution)