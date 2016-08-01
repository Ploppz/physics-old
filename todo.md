# Rendering ###
1. consistent names in Graphics and Renderers
2. FontRenderer should maybe render font independent of the scale (and rotation) of the scene

# Iterating ###
1. Polygon edge iterating: should have --. Maybe then it's more usable for general purpose.
  But it should also be random access. `set_to(int index)`. Definitely copy assignment operator. Not sure if operator= (int).
2. What if EdgeIterators contain some sort of EdgePointer or Edge, which has cached values of the end points, and index.
   It knows how to ++/--. But EdgeIterator takes care of "ending" and thus properly iterating.

... should be a goal that the `*iterator` should point to a struct _used by the iteration_ - that way it doesn't matter
    that there's some extra data


# Logical/Algorithm ###
0. CRITICAL: the closest edge of a vertex: occationally it won't move the polygons apart......
  - Let adjacent edges be "out directions". At a reference edge, the "out directions" must both point _out_ of the reference polygon.
0. semi-CRITICAL: what to do when are no vertices inside of each other?
## these two...:
 - I think 1st, it fails because all vertices are intersects - in this case, use LSS method.
 - Approx fails sometimes, so employ a new approx method:

1. Chain of collisions - Needed. 30 bodies inside bounding box with restitution = 0.5.
 * alternatively, we only need to decrease the number of collisions by finding a good order to resolve in
3. What to do when an intersection consists only of intersects?
4. What to do when an intersection has only one intersect that's like on the other side
5. Are bounding boxes updated mid-resolution / per iteration?
6. Tree iteration: maybe only continue on the branch if there was an actual intersection.
  - maybe this way, Breadth-first will work well.

7. HybridVertex 'edge1' and 'edge2' are analogous to polygons 'p' and 'q' sent as arguments to extract intersections.. could change names
    - maybe everything could be array with two elements??? that's way more dynamic. Furthermore: can we maybe always assume the right
       order of the two polygons? Depends on algorithms.
       - then we wouldn't need those "get_X_of_parent".. 

---
3. and 4. happen even when we have sufficient # iterations (but not as frequently)
4. is why untangling isn't perfect

# Programmatical ###
6. Intersection should have same iteration mechanism as Polygon



# Ideas ###

##### NOTEBOOK

Intersection resolution with Sides:
* a Side is contiguous part of the intersection which is part of only one polygonn.

#####

## Tree:
* HOW to iterate tree using PairOrderer::begin,end?
* in doubt whether Pair should not be `std::pair<int, int>`...

### Tree - Thoughts
* "displace distribution" is by mass - makes sense to sort by mass
  - but of course, velocity contributes to the depth
* when I think about sorting by relative velocity somehow... association: persistent contact... dunnnno

### Multiple intersections between two bodies
* Natural part of the tree.
* Can we solve it in a batch?

## temporal coherency in intersections
* is there a way we can more easily recalculate the intersections, when the polygons have only translated a bit since last time?

## optimize `linear_find_contact`.
 - search for sweepline intersection in the edges incident in the collision first


# Questions to self
* considering cache: would it be wise to actually keep all the Body data sorted according to which order we have updated the bodies?
 - would be a lot of (uncachefriendly) moving around..?
 - (sort seldomly)
## Nature of Force Propagation
* depending on order of resolution, it may happen this frame or the next - no big deal?












## (old problems?) ##

- Intersection extraction: When 2 Intersections have Intersects on the same edge, the algorithm only makes 1!


## Makefile ##
-
  rm .depend
  make
  .. why does it remake all targets?
- altering TimeResolutionAlg.h -> lots of recompilation
