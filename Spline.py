import numpy as np
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
#   3D plot of the
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, PchipInterpolator, interp1d
from math import floor, ceil

points = ((
    (0,0),
    (4.8, 14.3),
    (4.8*2.0, 17.7),
    (4.8*3.0, 19.2),
    (4.8*4.0, 20.0),
    (4.8*5.0, 20.3),
    (4.8*6.0, 20.1),
    (4.8*7.0, 19.8),
    (4.8*8.0, 19.1),
    (4.8*9.0, 18.4),
    (4.8*10.0, 17.2),
    (4.8*11.0, 16.0),
    (4.8*12.0, 14.7),
    (4.8*13.0, 13.1),
    (4.8*14.0, 10.8),
    (4.8*15.0, 7.7),
    (4.8*16.0, 0)
))
print(points)

spline = PchipInterpolator(*zip(*points))
spline2 = interp1d(*zip(*points), 'cubic')
x,y = zip(*points)
#Get points over spline
X_ = np.linspace(floor(min(x)), floor(max(x)), ceil(max(x)) * 10 - floor(min(x)) * 10)
Y_ = spline(X_)
Y_2 = spline2(X_)

print(list(zip(*points)))
# plot 
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(*zip(*points), 'ro')
ax.plot(X_, Y_)
ax.plot(X_, Y_2, 'r')
ax.set_xlim(-10, 80)
ax.set_ylim(-25, 25)
ax.set_xlabel('$x$ (cm)',fontsize=16)
ax.set_ylabel('\n$y$ (cm)',fontsize=16)
plt.show()