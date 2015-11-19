#pragma once

#include <vector>
#include <set>
#include <iostream>

#include <limits>
#include <climits>
#include <assert.h>
#include <utility>
#include <cmath>
#include <algorithm>

#include <GL/glew.h>

/** Some notes about the implementation:
		- Coordinate system has origin bottom-left (also the input rectangles)	
		- No rotation of rectangles
		- Keeps a simple datastructure: a skyline
		- First estimates the final width, and sets it. Height is varaible - the packaging only expands
			upwards, like tetris.
*/

namespace SBL {
	// Area
	template <typename T>
	struct Area {
		T key;
		double width, height, x, y;
	};

	// Size
	template <typename T>
	struct Size {
		Size(double width, double height, T key) { this->width = width; this->height = height; this->key = key;}
		Size() {};
		T key;
		double width, height;
	};

	// Point
	struct Point
	{
		double x, y;

		Point(double x, double y) {
			this->x = x;
			this->y = y;
		}
		Point() {}

		bool operator< (const Point &rhs) const {
			return (x < rhs.x);
		}
	};
	struct PointLessThan {
		bool operator() (const Point &lhs, const Point &rhs)
		{
			return lhs.x < rhs.x;
		}
	};

	typedef std::set<Point> Skyline;
	typedef std::set<Point>::iterator skyline_ptr;

	// OUTPUT

	template <typename T>
	std::ostream &operator<< (std::ostream &lhs, const Size<T> &rhs)
	{
		lhs << "Size(" << rhs.width << ", " << rhs.height << ")";
		return lhs;
	}
	template <typename T>
	std::ostream &operator<< (std::ostream &lhs, const Area<T> &rhs)
	{
		lhs << "Area(" << rhs.width << ", " << rhs.height << ") at (" << rhs.x << ", " << rhs.y << ")\t\t key: " << rhs.key;
		return lhs;
	}
	std::ostream &operator<< (std::ostream &lhs, const Point &rhs);
	template <typename T>
	std::ostream &operator<< (std::ostream &lhs, const std::vector<Area<T>> &rhs)
	{
		std::cout << "Atlas:" << std::endl;
		for (auto it = rhs.begin(); it != rhs.end(); it ++)
		{
			std::cout << "\t" << *it << std::endl;
		}
		return lhs;
	}

	std::ostream &operator<< (std::ostream &lhs, const std::set<Point> &rhs);

	template <typename T>
	bool rect_greaterThan(Size<T> rect1, Size<T> rect2) {
		return (rect1.height > rect2.height);
	}

	// Class

	template <typename T>
	class SkylineBL_alg
	{
	private:
		std::vector<Size<T>> rectangles;
		
		double width, height;
		Skyline skyline;

		/* ------------ */
		skyline_ptr findBestPosition(Size<T> rect)
		{
			double newHeight, minHeight = std::numeric_limits<double>::max();
			// Keep track of best position in terms of pointer into the skyline set.
			std::set<Point>::iterator best;
			bool bestFound = false;
			double bestHeight = INT_MAX;

			std::set<Point>::iterator forward;
			double currentHeight, currentWidth;
			for (skyline_ptr point = skyline.begin(); point != skyline.end(); point ++)
			{
				forward = point; forward ++;
				bool tooWide = (point->x + rect.width > width);
				if (tooWide) break;
				// Estimate new height
				currentHeight = totalHeightIfInserted(rect, point);
					
				if (currentHeight < bestHeight) {
					bestFound = true;
					best = point;
					bestHeight = currentHeight;
				}
			}
			assert(bestFound);

			return best;
		}



		/** Emulates insertion into the skyline at a given point */
		double totalHeightIfInserted(Size<T> rect, skyline_ptr where)
		{
			double shapeEnd = where->x + rect.width;
			skyline_ptr tallest = where;

			// Find tallest skyline point that is under this shape:
			for (skyline_ptr it = where; it != skyline.end(); it ++)
			{
				if (it->x >= shapeEnd) break;	
				if (it->y > tallest->y) tallest = it;
			}
			return tallest->y + rect.height;

		}
		/**	Inserts rect to our atlas, and updates skyline*/
		void insert(Size<T> rect, skyline_ptr where, std::vector<Area<T>> &atlas)
		{
			Point where_p = *where;

			// Project end down to get the new skyline point.
			Point end;
			end.x = where_p.x + rect.width;
			end.y = project(end.x, skyline.begin(), skyline.end());

			// Remove all points in the skyline, under this rectangle
				// At the same time, find the height in which we will place this rect. (same functionality as
				// totalHeightIfInserted)
			double removeUntil = where->x + rect.width;
			double maxHeight = 0;
			skyline_ptr next;
			for (skyline_ptr removal = where; removal != skyline.end(); )
			{
				next = removal; next ++;
				if (removal->x <= removeUntil && removal->x >= where_p.x) {
					// Calculating maximum height
					// @ The second condition: if removal->x == removeUntil, we don't want to count that point
					if (removal->y > maxHeight && removal->x < removeUntil) maxHeight = removal->y;
					// Remove.
					removal = skyline.erase(removal);
				} else {
					// Remove no more
					break;
				}
			}
			// Add two new points: At start (upper left) and at end (bottom right) of rectangle.
			// To find the end point, we need to project the end of the rectangle down to the skyline:

			
			Point skylinePoint1(where_p.x, maxHeight + rect.height);
			skyline.insert(skylinePoint1);
			Point skylinePoint2(end.x, end.y);
			skyline.insert(skylinePoint2);

			// Insert into atlas
			Area<T> newArea;
			newArea.x = where_p.x; newArea.y = maxHeight;
			newArea.width = rect.width; newArea.height = rect.height;
			newArea.key = rect.key;

			atlas.push_back(newArea);
			// Update width && height.

			height = std::max(height, skylinePoint1.y);
		}




		// This function finds the y-value of the skyline given an x-value
		double project(double x, skyline_ptr iterationBegin, skyline_ptr iterationEnd)
		{
			// Loop through skyline
			skyline_ptr next;
			for (skyline_ptr it = iterationBegin; it != iterationEnd; it ++)
			{
				next = it; next ++;
				if (next == iterationEnd) return it->y;
				// Check if `x` lies within this point and the next point
				if (x >= it->x && x < next->x)
				{
					return it->y;
				}
			}	
			// I REALLY DON'T KNOW WHAT TO DO HERE
			std::cerr << "Skyline_alg::project: Projection returned 0 because it really had no clue." << std::endl;
			return 0;
		}


		void sortByHeight(std::vector<Size<T>> &rects)
		{
			std::sort(rects.begin(), rects.end(), rect_greaterThan<T>);
		}

	public:
		SkylineBL_alg(std::vector<Size<T>> rectangles)
		{
			this->rectangles = rectangles;
			skyline = {};
			width = height = 0;
		}



		std::vector<Area<T>> compute(double &outWidth, double &outHeight)
		{
			std::vector<Area<T>> atlas = {};
			width = height = 0;

			// Estimate width:
			{
				double totalArea = 0;
				for (auto it = rectangles.begin(); it != rectangles.end(); it ++)
				{
					totalArea += it->width * it->height;
				}
				// sqrt(average area)
				width = sqrt(totalArea*1.1);
			}
			// Sort by height.
			sortByHeight(rectangles);

			// Insert first point:
			skyline.insert(Point(0, 0));

			// Loop through input rectangles
			for (auto rect = rectangles.begin(); rect != rectangles.end(); rect ++)
			{

				// Find best place to put it:
				//  1. Minimize additional total height
				//  2. If tie: pick the space with the smallest maximum side 
				skyline_ptr best = findBestPosition(*rect);
				
				// `best` now holds the position IN SKYLINE at which we should place this box
				insert(*rect, best, atlas);

			}
			outWidth = width;
			outHeight = height;
			return atlas;
		}

	};




	/*
	 *  GL FUNCTIONS TO HELP DRAW RESULTS
	 */
	void addPoint(GLfloat *data, int &i, double x, double y, double r, double g, double b);

	template <typename T>
	GLfloat *generateVertices(const std::vector<Area<T>> &atlas, int &length)
	{
		/* Format:
		 * <x> <y> <r> <g> <b>
		 */
		// Need 6 vertices for each square
		GLfloat *r = new GLfloat[atlas.size() * 6 * 5];
		int i = 0;
		auto it = atlas.begin();
		for ( ;it != atlas.end(); it ++ )
		{
			// Note: addPoint alters i

			addPoint(r, i, it->x, 				it->y,				0, 0, 0);
			addPoint(r, i, it->x + it->width, 	it->y + it->height, .2f, .2f, .2f);
			addPoint(r, i, it->x, 				it->y + it->height, .2f, 0.5f, .2f);

			addPoint(r, i, it->x, 				it->y, 				0, 0.5f, 0);
			addPoint(r, i, it->x + it->width, 	it->y, 				0, 0.5f, 0);
			addPoint(r, i, it->x + it->width,	it->y + it->height, 0, 0.3f, 0);

		}

		length = atlas.size() * 6 * 5;

		return r;
	}
}

