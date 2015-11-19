#include "skyline-bl.h"

#include <iostream>


namespace SBL {
	std::ostream &operator<< (std::ostream &lhs, const Point &rhs)
	{
		lhs << "Point(" << rhs.x << ", " << rhs.y << ")";
		return lhs;
	}

	std::ostream &operator<< (std::ostream &lhs, const std::set<Point> &rhs)
	{
		std::cout << "Skyline: " << std::endl;
		for (auto it = rhs.begin(); it != rhs.end(); it ++)
		{
			std::cout << "\t" << *it << std::endl;
		}
		return lhs;
	}	

	void addPoint(GLfloat *data, int &i, double x, double y, double r, double g, double b)
	{
		// Temporary solution: divide by 500.
		x = x / 50 - 0.5;
		y  = y / 50 - 0.5;
		data[i] = static_cast<GLfloat>(x); ++ i;
		data[i] = static_cast<GLfloat>(y); ++ i;
		data[i] = static_cast<GLfloat>(r); ++ i;
		data[i] = static_cast<GLfloat>(g); ++ i;
		data[i] = static_cast<GLfloat>(b); ++ i;

	}
}


