## Process
Given a table of offsets representing the Height Above Baseline (HAB) and Half-Breadth (HB) of the keel, chines, gunwale, and deckridge of a kayak at defined station distances along the length of the kayak. Treat the HAB as the Z dimension, HB as the X dimension, and the station distamce as the Y dimension.

Additional inputs are:
1. Frame material thickness
2. Keel frame width
3. Transverse frame width

### Define stringer shapes
#### Chines and gunwale
1. Find the plane that best fits the (X,Y,Z) 3d points for that chine or gunwale.
2. Project each of those points to the best-fit plane.
3. Calculate the distance between each original on projected point and ensure it is within a reasonable tolerance (1mm?)
4. Find the arc that best fits the points
5. find the points where that arc would intersect the YZ plane.These are the bow and stern endpoints for that chine/gunwale. Add them to the list of points.
6. Interpolate a minimum-energy b-spline through the points

#### Keel
1. Find  the arc that best fits the keel points
2. Find where that arc intersects the best-fit plane of the chine with the lowest Z values. These are the endpoints of the keel.
3. Add those endpoints to the list of keel points
4. Interpolate a minimum-energy b-spline through those points

#### Deckridge
1. Find the two nearly-linear segments starting from the bow and stern of the kayak
2. For each segment, extrapolate to where it has the same y-value as the first and/or last point of the gunwale
3. Add those points to the appropriate deckridge segment
4. Represent each deckridge segment as a straight line between the first and last points of that segment.

### Draw transverse frames
1. Define a plane parallel to the XZ plane at the station distance
2. Find the intersection of each chine, the gunwale, the keel, and the deckridge with that plane
3. Draw a horizontal line rightward starting at (0, keel HAB at station + keel_frame_width) with to (0.25 * material_thickness, keel HAB at station + keel frame width)
4. Draw a vertical line downward with length keel_frame_width - 0.25 * material_thickness
5. Construct a line indicating the intersection of chine_plane_n with the station plane. Construct a line that is 0.25 * material_thickness along this line toward the vertical axis, perpendicular to the chine plane line. Call this chine_n_inset_line.
6. Find the point that is 0.25 * material_thickness along chine_n_inset_line, toward the X axis.
7. Draw an arc from the endpoint of the line in step 4 to the point found in step 6
8. Draw a line parallel to the chine plane instersection line, with length chine_depth - 0.25 * material_thickness, toward the Y axis.
9. Draw a line parallel to chine_n_inset_line, with length material_thickness, away from the X and Y axes.
10. Draw a line parallel to chine plane intersection line, with length chine_depth - 0.25 * material_thickness, away from the Y axis.
11. Repeat 5 - 10 for each additional chine.
12. Perform 5 - 9 for the gunwale
13. Draw a line to (deckridge HB at station + 0.5 * material*thickness, deckridge HAB at station).
14. Draw a vertical line downward of length deckridge_frame_width - 0.25 * material_thickness
15. Draw a horizontal line to the left of length material_thickness. If the line intersects the Y axis, end the line at the Y axis.
16. If the line did not intersect the Y axis, draw a vertical line upward of length deckridge_material_frame_width - 0.25 * material_thickness
17. If currently at x > 0, draw a horizontal line to the left until intersecting the Y axis.
18. Mirror all of the lines drawn above across the Y axis.

### Draw keel/deckridge frame
1. Find the greater of the Z-value of the gunwale spline or the deckridge line segment at the bow Y-value (where the gunwale spline intersects the YZ plane)
2. Do the same for the stern
3. The outline of the keel/deckridge frame consists of the first deckridge segment, a line from the bow point found in step 1 to the endpoint of the keel spline, the keel spline, a line from the other endpoint of the keel spline to the stern endpont from step 2, the second deckridge line segment and the extrapolation of that line forward to the global y-value of the rearmost point of the first deckridge segment, then vertically to connect with the first deckridge segment.
4. For the inner outline of the frame, make an offset of the outer outline toward the inside. Offset by a configurable distance.
5. Draw a rectangular notch centered at each station location in the top edge of both the upper and lower parts of the frame. The notch is material_thickness + slot_tolerance wide and 1/2 frame_thickness - 0.25 * material_thickness deep

### Draw stringer shapes
1. For each chine and the gunwale, draw the b-spline for that item in the plane for that item. That is, draw it as a 2D curve in its plane.
2. Offset the b-spline by chine_thicknedd or gunwale_thickness toward the inside of the b-spline curve
3. find the intersection of the YZ plane with the chine/gunwale plane. Draw this line in the two locations where it connects the b-spline curve to the offset one.
4. Find the intersecton of each station plane with this chine/gunwale plane. Draw a rectangular notch centered on each intersection line, with width material_thickness + slot_tolerance, and depth 1/2 frame_thickness - 0.25 * material_thickness

### Outputs to save
1. B-splines of each chine, the gunwale, keel. Plane (point + normal) for each chine and gunwale. Deckridge line segments. (Format TBD)
2. SVG for each station
3. SVG for keel/deckridge unit, each chine, and gunwale