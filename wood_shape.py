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
    
    def calculate_strain_energy(self):
        """Calculate approximate strain energy in the bent wood"""
        t = np.linspace(0, 1, self.resolution)
        curve = self.spline(t)
        # Second derivative approximates curvature for small deflections
        second_deriv = self.spline.derivative(2)(t)
        return np.sum(np.square(second_deriv)) * self.tension
    
    def draw(self, filename, width=800, height=600, margin=50):
        curve_points = self.get_curve_points()
        
        # Scale points to fit in the drawing area
        x_min, y_min = np.min(curve_points, axis=0)
        x_max, y_max = np.max(curve_points, axis=0)
        
        scale_x = (width - 2 * margin) / (x_max - x_min)
        scale_y = (height - 2 * margin) / (y_max - y_min)
        scale = min(scale_x, scale_y)
        
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
        ctx = cairo.Context(surface)
        
        # Fill background
        ctx.set_source_rgb(1, 1, 1)
        ctx.paint()
        
        # Draw the curve
        ctx.set_source_rgb(0.4, 0.2, 0)  # Wood-like brown color
        ctx.set_line_width(3)
        
        # Transform coordinates to fit the drawing area
        for i, point in enumerate(curve_points):
            x = margin + (point[0] - x_min) * scale
            y = margin + (point[1] - y_min) * scale
            if i == 0:
                ctx.move_to(x, y)
            else:
                ctx.line_to(x, y)
        
        ctx.stroke()
        
        # Draw control points
        ctx.set_source_rgb(0.8, 0, 0)
        for point in self.points:
            x = margin + (point[0] - x_min) * scale
            y = margin + (point[1] - y_min) * scale
            ctx.arc(x, y, 4, 0, 2 * np.pi)
            ctx.fill()
        
        surface.write_to_png(filename)

def generate_bent_wood_shape(points, tension=1.0, resolution=100):
    """Create a bent wood shape from a list of points"""
    return BentWoodShape(points, tension, resolution)