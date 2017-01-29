/*
 * test_semi_plan_wall.cpp
 *
 *  Created on: Jan 28, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#include "voro++.hh"

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <chrono>
#include <ctime>
#include <unordered_map>
#include <cmath>

#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/lagged_fibonacci.hpp>

#include <numeric>

struct wall_half_plane : public voro::wall {
public:
	/** Constructs a half-plane wall object.
	 * \param[in] (xc_,yc_,zc_) a normal vector to the plane.
	 * \param[in] ac_ a displacement along the normal vector.
	 * \param[in] w_id_ an ID number to associate with the wall for
	 *                  neighbor tracking.
	 * For now, it only does half-planes with a cutoff defined on the
	 * X direction...
	 *
	 */


	wall_half_plane(double xc_,double yc_,double zc_,double ac_, double cutoffx_, int cutoffxsign_,  double cutoffz_, int cutoffzsign_, int w_id_=-99)
					: w_id(w_id_), xc(xc_), yc(yc_), zc(zc_), ac(ac_),
					  cutoffx(cutoffx_), cutoffxsign(cutoffxsign_),
					  cutoffz(cutoffz_), cutoffzsign(cutoffzsign_) {}

	bool point_inside(double x,double y,double z)
	{
		return (x*xc+y*yc+z*zc<ac) && (x>cutoffx*cutoffxsign);
	};

	template<class v_cell>
	bool cut_cell_base(v_cell &c,double x,double y,double z)
	{
		double tol = 0.01;
		double distancex = x - cutoffx + tol * cutoffxsign;
		if(distancex * cutoffxsign > 0)
		{
	    	double dq=2*(ac-x*xc-y*yc-z*zc);
	    	return c.nplane(xc,yc,zc,dq,w_id);
		}
	    return true;
	}

	bool cut_cell(voro::voronoicell &c,double x,double y,double z) {return cut_cell_base(c,x,y,z);}
   	bool cut_cell(voro::voronoicell_neighbor &c,double x,double y,double z) {return cut_cell_base(c,x,y,z);}

private:

	const int w_id;
	const double xc,yc,zc,ac;
	const double cutoffx;
	const int cutoffxsign;
	const double cutoffz;
	const int cutoffzsign;
};

int main()
{
	double x_min,x_max,y_min,y_max,z_min,z_max;
	x_min = -1;
	x_max = 1;
	y_min = -0.5;
	y_max = 0.5;
	z_min = -0.5;
	z_max = 0.5;

	voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,20,10,10,
                        false,false,false,8);

	wall_half_plane upper_cut(0.1,0,-1,0,
								0,1,0,1,10);
	wall_half_plane lower_cut(0.1,0,1,0,
								0,1,0,-1,10);

	con.add_wall(upper_cut);
//	con.add_wall(lower_cut);

	boost::random::lagged_fibonacci607 m_rng;
	boost::random::uniform_01<> rnd;

	double x,y,z;
	for(int iii=0; iii<20; iii++) {
		x=x_min+rnd(m_rng)*(x_max-x_min);
		y=y_min+rnd(m_rng)*(y_max-y_min);
		z=z_min+rnd(m_rng)*(z_max-z_min);
		con.put(iii,x,y,z);
	}

	con.draw_cells_gnuplot("test_both_upper.dat");
}
