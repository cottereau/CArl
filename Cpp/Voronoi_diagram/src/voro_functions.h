/*
 * voro_functions.h
 *
 *  Created on: Jan 29, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef VORO_FUNCTIONS_H_
#define VORO_FUNCTIONS_H_

// Voronoi calculation example code
#include "voro++.hh"
#include "getpot.h"
#include <stdlib.h>

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

#define BOUNDARY_ID_MIN_Z 1
#define BOUNDARY_ID_MIN_Y 2
#define BOUNDARY_ID_MAX_X 3
#define BOUNDARY_ID_MAX_Y 4
#define BOUNDARY_ID_MIN_X 5
#define BOUNDARY_ID_MAX_Z 6

// Point hash function
struct PointHash_3D {
	std::size_t operator()(const std::vector<long>& k) const
	{
		long prime0 = 73856093;
		long prime1 = 19349669;
		long prime2 = 83492791;
		long primeN = 2038074743;

		return ( ( k[0] * prime0 ) ^ ( k[1] * prime1 ) ^ ( k[2] * prime2 ) ) % primeN;
	}
};

struct PointHash_3D_Equal {
	bool operator()(const std::vector<long>& lhs, const std::vector<long>& rhs) const
	{
		return lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2];
	}
};

class point_hasher
{
private:
	std::vector<double> m_Grid_MinPoint;
	std::vector<double> m_Grid_MaxPoint;

	double m_eps;
	long m_GridN_min;
public:

	point_hasher(long GridN_min = static_cast<long>(1E9)) : m_GridN_min (GridN_min)
	{
		m_eps = -1;
		m_Grid_MinPoint.resize(3);
		m_Grid_MaxPoint.resize(3);
	};

	void set_grid_constraints(const std::vector<double>& minPoint, const std::vector<double>& maxPoint)
	{
		m_Grid_MinPoint = minPoint;
		m_Grid_MaxPoint = maxPoint;

		std::vector<double> eps_candidates(3,0);
		for(unsigned int iii = 0; iii < 3; ++iii)
		{
			eps_candidates[iii] = (m_Grid_MaxPoint[iii] - m_Grid_MinPoint[iii])/m_GridN_min;
		}

		m_eps = *std::min_element(eps_candidates.begin(),eps_candidates.end());

		for(unsigned int iii = 0; iii < 3; ++iii)
		{
			m_Grid_MinPoint[iii] -= 2*m_eps;
			m_Grid_MaxPoint[iii] += 2*m_eps;
		}
	}

	void convert_to_discrete(const std::vector<double>& iPoint, std::vector<long>& oPoint)
	{
		oPoint[0] = lround( (iPoint[0] -  m_Grid_MinPoint[0] )/m_eps);
		oPoint[1] = lround( (iPoint[1] -  m_Grid_MinPoint[1] )/m_eps);
		oPoint[2] = lround( (iPoint[2] -  m_Grid_MinPoint[2] )/m_eps);
	}
};

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


	wall_half_plane(double xc_,double yc_,double zc_,double ac_, double cutoffx_, int cutoffxsign_,  int w_id_=-99)
					: w_id(w_id_), xc(xc_), yc(yc_), zc(zc_), ac(ac_),
					  cutoffx(cutoffx_), cutoffxsign(cutoffxsign_) {}

	bool point_inside(double x,double y,double z)
	{
		return (x*xc+y*yc+z*zc<ac) && (x>cutoffx*cutoffxsign);
	};

	template<class v_cell>
	bool cut_cell_base(v_cell &c,double x,double y,double z)
	{
		double tol = 0.01;
		double distancex = x - cutoffx + tol * cutoffxsign;
		if(distancex * cutoffxsign > 0 )

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
};

bool IsNearBorder(double px, double py, double pz, double a, double b, double c, double d, double eps);

void voro_statistics(voro::container &con, double& nb_of_faces, double& nb_of_edges);

void NewExportVoronoiToGmsh(voro::container &con, double weight, std::string &baseFilename);

void ExportVoronoiToGmshTest(voro::container &con, std::string &outFilename);

#endif /* VORO_FUNCTIONS_H_ */
