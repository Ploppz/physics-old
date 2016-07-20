# Rendering
1. consistent names in Graphics and Renderers
2. FontRenderer should maybe render font independent of the scale (and rotation) of the scene

## old thoughts about it
* Separate GL interfacing (uploading, rendering) and the data construction
* Persistent buffers of e.g. glyph vertex data
* FontTexture should upload the data and return the GL ID
* some Buffer class for `pos 2 tex 2 col 3` - includes texture
  - no, we need a buffer specialized for Font Rendering, which contains a FontTexture
  --> Should all these different specializations inherit from something common? maybe just for polymorphism - BufferObject - VisualData
    - this way, we can put all Buffer Objects into one std::vector

* TAKE1 Renderer has all interfaces, and contains the data. While other object such as LineBuffer/FontRenderer, know
    how to fill those buffers, and how to upload/render them
  - compound draw operations = draw operations to single interfaces (LineBuffer, FontRenderer)
  * OR:
    - Draw interfaces have buffers, and `get_buffer(std::string)`
    - .. FontRenderer needs to have `add_text`
* TAKE2
  - Renderer is an abstract class 
    - `set_active_buffers(std::string)`
    - contains buffers, knows how to add to them and how to render them
  - Graphics: Grand interface which contains (explicitly listed) Renderer instances
  - `Graphics::set_active_buffers(std::string)`

* IDEA: FontRenderer and LineRenderer can stay independent of e.g. Polygon:
  -  LineRenderer only draws lines.. Graphics tells LineRenderer which lines to draw.

# Logical/System
0. Sort the indices of one single pair
1. Chain of collisions - Needed. 30 bodies inside bounding box with restitution = 0.5.
 * alternatively, we only need to decrease the number of collisions by finding a good order to resolve in
2. Depth distribution - probably okay now with only mass
3. What to do when an intersection consists only of intersects?
4. What to do when an intersection has only one intersect that's like on the other side
5. Are bounding boxes updated mid-resolution / per iteration?

---
3. and 4. happen even when we have sufficient # iterations (but not as frequently)
4. is why untangling isn't perfect

# Programmatical
6. Intersection should have same iteration mechanism as Polygon



# Ideas

## Intersection Tree
* may tell which pair to resolve first - based on mass or intersection depth or velocity - a *mix* of these
 - maybe sum(mass)
* in the case of intersection depth, we may reuse some of the intersection calculations used
* maybe it could be good to just make sure to start somewhere else every time?

### Multiple intersections between two bodies
* Natural part of the tree.
* Can we solve it in a batch?

### thought
* is there a way we can more easily recalculate the intersections, when the polygons have only translated a bit since last time?

## Intersection Tree 2
* create a tree from the pairs of SAP.
  - each Pair is one node
  - create edges between Pairs that share one box-id
    - hashmap [box-id] -> std::vector<Pair>
* sort edges of the tree by sum of mass
  - O
* traverse (will need some 'processed' bool) breadth-first --- and treat along edges

* goal:
  - take pair with highest mass first, then continue from there..
### IN PROGRESS:
  - HOW to iterate tree using PairOrderer::begin,end?
  - in doubt whether Pair should not be `std::pair<int, int>`...



### Intersection Tree 2 - Thoughts
 - "displace distribution" is by mass - makes sense to sort by mass
    - but of course, velocity contributes to the depth
 - when I think about sorting by relative velocity somehow... association: persistent contact... dunnnno




## More performant intersection extraction

{
    - Instead of int.ext., we may intersect polygon A with sweep lines of polygon B (and vice versa...)?
    - pros: continuous.
    - cons: quite expensive for an overlap test.
    - cons: can't group vertices in intersections
}
* Instead, we may try to optimize `linear_find_contact`.
 - search for sweepline intersection in the edges incident in the collision first


# Questions to self
* considering cache: would it be wise to actually keep all the Body data sorted according to which order we have updated the bodies?
 - would be a lot of (uncachefriendly) moving around..?
## Nature of Force Propagation
* depending on order of resolution, it may happen this frame or the next - no big deal?












###### lalala #########################

- Crashes when one polygon is OUTSIDE mode
  - before anything else:`get_contact` needs to properly decide which direction each LineStripSeries::Vertex should go
- Intersection extraction: When 2 Intersections have Intersects on the same edge, the algorithm only makes 1!


- 4 Intersects:
  - not sure but I think shadow casting is off  --- might appear like it casts the wrong way?
  - ... OF COURSE: normal is still way wrong


## Makefile ##
-
  rm .depend
  make
  .. why does it remake all targets?
- altering TimeResolutionAlg.h -> lots of recompilation
