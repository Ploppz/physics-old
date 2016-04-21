# Rendering Polygons #
 - Indices rather than vertices


# Iterating a Polygon #

# Logic #
- `extract_intersections` takes flip_q and flip_p (rename to insideout or something?)
    - `calculate_first_contact` (or rather `rewind_out_of`) should also flip the inside test
    - We need some nice system for this.. e.g.  each polygon has a _level_
- ABSOLUTE position type:
    - cannot change position: it acts like a hole in containing polygon
    - Impulses just propagate to parent polygon, if any

# Other #
Suspicion: alpha in how-near-edge


# Collision resolution discrepancies #
 - if there was a collision last frame, 

 ... does the physical response guarantee that the points will be separating?

