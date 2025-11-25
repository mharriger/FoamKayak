from stations_loader import load_stations_file
import skspatial.objects as skso
from minimum_energy_bspline import minimum_energy_bspline
import numpy as np

from freecad_functions import (makeFreeCADDocument, add_bspline_sketch)

from geom_functions import (
    global_to_local,
    local_to_global,
    make_coordinate_system,
    printPointDistances,
    find_chine_endpoints,
    segment_polyline_near_straight,
    fit_plane,
    project_points_to_plane,
    planarize_and_extrapolate_chine,
    offset_point_along_plane_intersection
    )



#TODO: Better b-spline representation (knot vector, multiplicities, control points)

KAYAK_NAME = 'SeaRoverST'


data = load_stations_file(f'data/{KAYAK_NAME}.offsets')
stations = data['stations']
keelZ = data['keel_hab']
chines = data['chines']
if len(chines) == 0:
    chinesX = np.empty((0, len(stations)))
    chinesZ = np.empty((0, len(stations)))
else:
    chinesX = np.vstack([c['hb'] for c in chines])
    chinesZ = np.vstack([c['hab'] for c in chines])
deckridgeZ = data['deckridge']['hab']
gunwaleX = data['gunwale']['hb']
gunwaleZ = data['gunwale']['hab']

if chinesX.ndim == 1:
    chinesX.shape = (1, len(chinesX))
if chinesZ.ndim == 1:
    chinesZ.shape = (1, len(chinesZ))

chine_bsplines = []
chine_planes = []
chine_points_list_2d = []
plane_coords = []

geom_dict = {}

def process_chine(chineX, chineZ):
    points = skso.Points([[chineX[i], stations[i], chineZ[i]] for i in range(len(stations))])

    # Make the chine points coplanar, and approximate the bow and stern endpoints
    threedpoints, twodpoints, chine_plane, local_coords, distances = planarize_and_extrapolate_chine(points)

    #TODO: Check that point distances are acceptable
    chine_bspline = minimum_energy_bspline(twodpoints)
    return threedpoints, twodpoints, chine_bspline, chine_plane, local_coords, distances

for idx in range(len(chinesX)):
    chineX = chinesX[idx]
    chineZ = chinesZ[idx]

    chine_points_3d, chine_points_2d, chine_bspline, chine_plane, local_coords, distances = process_chine(chineX, chineZ)
    chine_points_list_2d.append(chine_points_2d)
    chine_bsplines.append(chine_bspline)
    chine_planes.append(chine_plane)
    plane_coords.append(local_coords)
    geom_dict[f'Chine{idx+1}'] = {
        '3d_points': chine_points_3d,
        '2d_points': chine_points_2d,
        'bspline': chine_bspline,
        'plane_normal': chine_plane.normal,
        'plane_point': local_coords[0],
        'distances': distances
    }

gunwale_points_3d, gunwale_points_2d, gunwale_bspline, gunwale_plane, gunwale_coords, gunwale_distances = process_chine(gunwaleX, gunwaleZ)
geom_dict['Gunwale'] = {
    '3d_points': gunwale_points_3d,
    '2d_points': gunwale_points_2d,
    'bspline': gunwale_bspline,
    'plane_normal': gunwale_plane.normal,
    'plane_point': gunwale_coords[0],
    'distances': gunwale_distances
}


geom_dict['Keel'] = {
    '3d_points': [skso.Point([0, stations[i], keelZ[i]]) for i in range(len(stations))],
    '2d_points': [[stations[i], keelZ[i]] for i in range(len(stations))],
    'bspline': minimum_energy_bspline([[stations[i], keelZ[i]] for i in range(len(stations))]),
    'plane_normal': np.array([1,0,0]),
    'plane_point': skso.Point([0,0,0]),
    'distances': [0 for _ in range(len(stations))]
}

doc = makeFreeCADDocument(KAYAK_NAME)

for idx, chine_bspline in enumerate(chine_bsplines):
    chine_plane = chine_planes[idx]
    chine_coords = plane_coords[idx]
    add_bspline_sketch(doc, f'Chine{idx+1}', chine_coords[0], chine_plane.normal, chine_bspline, chine_points_list_2d[idx])
add_bspline_sketch(doc, 'Gunwale', gunwale_coords[0], gunwale_plane.normal, gunwale_bspline, gunwale_points_2d)

add_bspline_sketch(doc, 'Keel', skso.Point([0,0,0]), np.array([1,0,0]), geom_dict['Keel']['bspline'], geom_dict['Keel']['2d_points'], rotation=(0,1,0,270))

doc.recompute()
doc.saveAs(f'{KAYAK_NAME}_bsplines.FCStd')