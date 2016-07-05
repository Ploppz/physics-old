# Currently WIP #

1. Depth distribution
2. Chain of collisions

- Consider making a Placement struct instead of PlacementBackup, which also e.g. Polygon can use, and Body::placement() can return




- Crashes when one polygon is OUTSIDE mode
  - before anything else: get_contact needs to properly decide which direction each LineStripSeries::Vertex should go
- Intersection extraction: When 2 Intersections have Intersects on the same edge, the algorithm only makes 1!


- 4 Intersects:
  - not sure but I think shadow casting is off  --- might appear like it casts the wrong way?
  - ... OF COURSE: normal is still way wrong


  .. MakeFile:
  rm .depend
  make
  .. why does it remake all targets?
  - altering TimeResolutionAlg.h -> lots of recompilation

# Rendering Polygons #
 - Indices rather than vertices

# Polygon calculations (area, etc) #
ensure either that it updates variables automatically, or that it crashes 
- getters and setters

# Iterating a Polygon #

# Logic #
- `extract_intersections` takes `flip_q` and `flip_p` (rename to insideout or something?)
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

