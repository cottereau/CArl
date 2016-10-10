/*
 * search_point_2.h
 *
 *  Created on: Aug 25, 2015
 *      Author: Thiago Milanetto Schlittler
 *
 *     		Special point class, needed to save the vertex handles while doing
 *     	the search.
 */

#ifndef INTERSECTIONS_SEARCH_POINT_2_H_
#define INTERSECTIONS_SEARCH_POINT_2_H_

#include "carl_headers.h"
#include "triangular_mesh_2.h"

/*
 *  --- TreePoint class
 *
 *  	This is essentially a Point_2 class, but with an extra information, the
 *  "vHandle", a vertex handle used to connect it with the triangulation data
 *  structure. Based on CGAL's file "[CGAL root]/Spatial_searching/Point.h"
 *
 */
class TreePoint_2 {
	double vec[2];
	Vertex_handle_2	vHandle;

	TreePoint_2() { vec[0]= vec[1] = 0; vHandle = Vertex_handle_2();}
	TreePoint_2(double x, double y, Vertex_handle_2 inputHandle)
	{
		vec[0]=x;
		vec[1]=y;

		vHandle = inputHandle;
	}

	TreePoint_2(Vertex_handle_2 inputHandle)
	{
		vec[0]=inputHandle->point().x();
		vec[1]=inputHandle->point().y();

		vHandle = inputHandle;
	}

	double x() const { return vec[ 0 ]; }
	double y() const { return vec[ 1 ]; }
	Vertex_handle_2 handle() const { return vHandle; }

	double& x() { return vec[ 0 ]; }
	double& y() { return vec[ 1 ]; }

	void	set(Vertex_handle_2 inputHandle)
	{
		vec[0]=inputHandle->point().x();
		vec[1]=inputHandle->point().y();

		vHandle = inputHandle;
	}

	bool operator==(const TreePoint_2& p) const
	{
		return (x() == p.x()) && (y() == p.y()) ;
	}

	bool  operator!=(const TreePoint_2& p) const { return ! (*this == p); }

}; //end of class

namespace CGAL {
	  template <>
	  struct Kernel_traits<TreePoint_2> {
			struct Kernel {
			  typedef double FT;
			  typedef double RT;
			};
	  };
}
struct Construct_coord_iterator {
	  typedef  const double* result_type;
	  const double* operator()(const TreePoint_2& p) const
	  { return static_cast<const double*>(p.vec); }
	  const double* operator()(const TreePoint_2& p, int)  const
	  { return static_cast<const double*>(p.vec+2); }
};



#endif /* INTERSECTIONS_SEARCH_POINT_2_H_ */
