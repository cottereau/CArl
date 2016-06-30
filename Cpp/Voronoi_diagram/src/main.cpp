//============================================================================
// Name        : Voronoi_tesselation.cpp
// Author      : Thiago Milanetto Schlittler
// Version     :
// Copyright   : Your copyright notice
// Description : Generation of a 3D polycristal structure
//				 using Voronoi diagrams
//============================================================================

// Voronoi calculation example code
#include "voro++.hh"
#include "getpot.h"

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

#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/lagged_fibonacci.hpp>

#define BOUNDARY_ID_MIN_Z 1
#define BOUNDARY_ID_MIN_Y 2
#define BOUNDARY_ID_MAX_X 3
#define BOUNDARY_ID_MAX_Y 4
#define BOUNDARY_ID_MIN_X 5
#define BOUNDARY_ID_MAX_Z 6

bool IsNearBorder(double px, double py, double pz, double a, double b, double c, double d, double eps)
{
	double difference = std::abs(px*a + py*b + pz*c - std::abs(d));

	if(difference < eps)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void ExportVoronoiToGmsh(voro::container &con, double weight, std::string &baseFilename)
{
	// ---------------
	//  Preamble
	// ---------------

	// Vector containing the vertices of each face
	/*
	 * --> The information of a face "k" with "n_k" vertices is located in the
	 *     middle of "face_vertices", with the structure
	 *
	 *     ... n_k, v_{0,n_k}, v_{1,n_k}, ..., v_{n_k - 1,n_k}, ...
	 *
	 *     where v_{i,n_k} is the i'th vertex. The positions of the "n_k"'s
	 *     follow a recursion
	 *
	 *     pos_{n_0} 		= 0
	 *     pos_{n_1} 		= (n_0 + 1)				= pos_{n_0} + (n_0 + 1)
	 *     pos_{n_2} 		= (n_0 + 1) + (n_1 + 1)	= pos_{n_1} + (n_1 + 1)
	 *
	 *     pos_{n_(k+1)}	= 						= pos_{n_k} + (n_k + 1)
	 *
	 *     we will use the auxiliary index "aux_index" to save it
	 */
	std::vector<int>	face_vertices_index;
	int aux_index;
	int face_nb_of_vertices;

	// Cell position
	double pos_x, pos_y, pos_z;

	// Loop variable used to run over all cells
	voro::c_loop_all	cell_loop(con);

	// Dummy cell with neighbour informations
	voro::voronoicell dummy_cell;

	// Some queues to save the boundaries
	std::queue<int> minXQueue;
	std::queue<int> maxXQueue;

	std::queue<int> minYQueue;
	std::queue<int> maxYQueue;

	std::queue<int> minZQueue;
	std::queue<int> maxZQueue;

	// -----------------------
	//  Gmsh data structure
	// -----------------------

	double eps = 1E-8*pow((con.bx-con.ax)*(con.by-con.ay)*(con.bz-con.az)/con.total_particles(),1./3.);
	std::cout << eps << std::endl;

	std::string  outFilename = baseFilename + ".geo";
	std::string  outFilenameCenters = baseFilename + "_centers.dat";

	std::ofstream gmshOutput(outFilename, std::ofstream::trunc);
	std::ofstream centerOutput(outFilenameCenters, std::ofstream::trunc);

	std::chrono::time_point<std::chrono::system_clock> time_now;
	time_now = std::chrono::system_clock::now();
	std::time_t timestamp = std::chrono::system_clock::to_time_t(time_now);

	gmshOutput << "// Gmsh file created using the Voronoi tesselation generator"
			   << std::endl
			   << "// Timestamp: " << std::ctime(&timestamp) << std::endl;

	gmshOutput << "// Number of cells : " << con.total_particles() << std::endl
			   << "// Dimensions      : (" << con.ax << ", " << con.bx << ") ("
			                               << con.ay << ", " << con.by << ") ("
										   << con.az << ", " << con.bz << ")"
			   << std::endl
	           << "// Volume          : " << con.sum_cell_volumes()
	           << std::endl << std::endl;

	gmshOutput << "// Weights" << std::endl;
    gmshOutput << "w1 = " << weight << ";" << std::endl << std::endl;

	// Vertices !
	std::vector<double> cell_vertices;
	int cell_nb_of_vertices = 0;

	// Index containing the global index of each vertex
	int vertex_global_index = 1;

	// Dummy boolean used as a marker to when to stop iterating
	bool keepIterating = cell_loop.start();

	while(keepIterating)
	{
		// Transfer needed data to "dummy_cell"
		con.compute_cell(dummy_cell,cell_loop);

		// Coordinates infos are saved in the local coord. system, must get the
		// 		cell particle's coords.

		cell_loop.pos(pos_x,pos_y,pos_z);

		dummy_cell.vertices(pos_x,pos_y,pos_z,cell_vertices);
		cell_nb_of_vertices = dummy_cell.p;

		for(int iii = 0; iii < cell_nb_of_vertices; ++iii)
		{
			gmshOutput 	<< "Point(" << vertex_global_index << ") = {"
						<< cell_vertices[3*iii] << ", "
						<< cell_vertices[3*iii+1] << ", "
						<< cell_vertices[3*iii+2] << ", w1};" << std::endl;
			++vertex_global_index;
		}

		// Go to next cell (returns "false" if it's over)
		keepIterating = cell_loop.inc();
	}

	gmshOutput << std::endl;

	// Edges and faces!
	int treated_vertices = 1;
	int treated_edges    = 1;
	int treated_faces    = 1;

	int cell_number_of_faces = 0;

	int current_vertex	 = 0;
	int current_edge     = 0;
	int current_face     = 0;

	int edge_start = 0;
	int edge_end = 0;

	int edge_global_index = 1;
	int face_global_index = 1;
	int cell_global_index = 1;

	double pointX = -1;
	double pointY = -1;
	double pointZ = -1;

	int minXCounter = 0;
	int maxXCounter = 0;

	int minYCounter = 0;
	int maxYCounter = 0;

	int minZCounter = 0;
	int maxZCounter = 0;

	keepIterating = cell_loop.start();

	centerOutput << con.total_particles() << std::endl;

	while(keepIterating)
	{
		// Transfer needed data to "dummy_cell"
		con.compute_cell(dummy_cell,cell_loop);

		// ... and extract it
		cell_loop.pos(pos_x,pos_y,pos_z);
		dummy_cell.vertices(pos_x,pos_y,pos_z,cell_vertices);

		// Get neighbors and face vertices
		dummy_cell.face_vertices(face_vertices_index);

		cell_number_of_faces = dummy_cell.number_of_faces();

		aux_index = 0;

		// Loop over all the faces
		for(int iii = 0; iii < cell_number_of_faces; ++iii)
		{
			face_nb_of_vertices = face_vertices_index[aux_index];


			// -> Extract the edges
			for(int jjj = 1; jjj < face_nb_of_vertices; ++jjj)
			{
				edge_start = treated_vertices + face_vertices_index[aux_index + jjj];
				edge_end   = treated_vertices + face_vertices_index[aux_index + jjj + 1];

				gmshOutput	<< "Line(" << edge_global_index << ") = {"
							<< edge_start << ", " << edge_end << "};"
							<< std::endl;
				++edge_global_index;
			}

			edge_start = treated_vertices + face_vertices_index[aux_index + face_nb_of_vertices];
			edge_end   = treated_vertices + face_vertices_index[aux_index + 1];

			gmshOutput	<< "Line(" << edge_global_index << ") = {"
						<< edge_start << ", " << edge_end << "};"
						<< std::endl;
			++edge_global_index;

			// -> Test if the face is on a border
			minXCounter = 0;
			maxXCounter = 0;

			minYCounter = 0;
			maxYCounter = 0;

			minZCounter = 0;
			maxZCounter = 0;

			for(int jjj = 1; jjj < face_nb_of_vertices + 1; ++jjj)
			{
				current_vertex = face_vertices_index[aux_index + jjj];

				pointX = cell_vertices[3*current_vertex];
				pointY = cell_vertices[3*current_vertex + 1];
				pointZ = cell_vertices[3*current_vertex + 2];

				// Test min X
				if(IsNearBorder(	pointX,pointY,pointZ,
									-1,0,0, con.ax, eps))
				{
					// Add to border min x
					++minXCounter;
				}

				// Test max X
				if(IsNearBorder(	pointX,pointY,pointZ,
									1,0,0, con.bx, eps))
				{
					// Add to border max x
					++maxXCounter;
				}

				// Test min Y
				if(IsNearBorder(	pointX,pointY,pointZ,
									0,-1,0, con.ay, eps))
				{
					// Add to border min y
					++minYCounter;
				}

				// Test max Y
				if(IsNearBorder(	pointX,pointY,pointZ,
									0,1,0, con.by, eps))
				{
					// Add to border max y
					++maxYCounter;
				}

				// Test min Z
				if(IsNearBorder(	pointX,pointY,pointZ,
									0,0,-1, con.az, eps))
				{
					// Add to border min z
					++minZCounter;
				}

				// Test max Z
				if(IsNearBorder(	pointX,pointY,pointZ,
									0,0,1, con.bz, eps))
				{
					// Add to border max z
					++maxZCounter;
				}
			}

			if(minXCounter == face_nb_of_vertices)
			{
				// It's on the min X face, add it!
				minXQueue.push(face_global_index);
			}

			if(maxXCounter == face_nb_of_vertices)
			{
				// It's on the max X face, add it!
				maxXQueue.push(face_global_index);
			}

			if(minYCounter == face_nb_of_vertices)
			{
				// It's on the min X face, add it!
				minYQueue.push(face_global_index);
			}

			if(maxYCounter == face_nb_of_vertices)
			{
				// It's on the max X face, add it!
				maxYQueue.push(face_global_index);
			}

			if(minZCounter == face_nb_of_vertices)
			{
				// It's on the min X face, add it!
				minZQueue.push(face_global_index);
			}

			if(maxZCounter == face_nb_of_vertices)
			{
				// It's on the max X face, add it!
				maxZQueue.push(face_global_index);
			}

			// -> Build the face
			gmshOutput	<< std::endl;
			gmshOutput	<< "Line Loop(" << face_global_index << ") = {";

			for(int jjj = 0; jjj < face_nb_of_vertices - 1; ++jjj)
			{
				current_edge = treated_edges + jjj;
				gmshOutput	<< current_edge << ", ";
			}

			// -> Close the face
			current_edge = treated_edges + face_nb_of_vertices - 1;
			gmshOutput	<< current_edge << "};" << std::endl;

			// -> Build the plane surface
			gmshOutput	<< "Plane Surface(" <<  face_global_index
						<< ") = {" <<  face_global_index << "};" << std::endl;
			gmshOutput	<< std::endl;

			// Update face index
			++face_global_index;
			treated_edges = edge_global_index;

			// -> Jump to the next face
			aux_index += face_nb_of_vertices + 1;
		}

		// -> Finally, build the cell!
		gmshOutput << "Surface Loop(" << cell_global_index << ") = {";
		for(int iii = 0; iii < cell_number_of_faces - 1; ++iii)
		{
			current_face = treated_faces + iii;
			gmshOutput	<< current_face << ", ";
		}

		// -> Close the cell
		current_face = treated_faces + cell_number_of_faces - 1;
		gmshOutput	<< current_face << "};" << std::endl;

		// -> Set up the volumes
		gmshOutput	<< "Volume(" <<  cell_global_index
					<< ") = {" <<  cell_global_index << "};" << std::endl;
		gmshOutput	<< "Physical Volume(" <<  cell_global_index
					<< ") = {" <<  cell_global_index << "};" << std::endl;
		gmshOutput	<< std::endl;

		// -> Print the center
		centerOutput 	<< pos_x << " " << pos_y << " " << pos_z << " "
						<< cell_global_index << std::endl;

		// Update the cell index
		++cell_global_index;
		treated_faces = face_global_index;

		// Go to next cell (returns "false" if it's over)
		treated_vertices += dummy_cell.p;
		keepIterating = cell_loop.inc();
	}

	gmshOutput << "Coherence;" << std::endl;

	// Now save the physical surfaces
	int dummySurface = -1;

	std::cout 		<< minXQueue.size() << " " << maxXQueue.size() << " "
					<< minYQueue.size() << " " << maxYQueue.size() << " "
					<< minZQueue.size() << " " << maxZQueue.size() << std::endl;

	gmshOutput 	<< "Physical Surface(" << BOUNDARY_ID_MIN_X << ") = {";
	while(!minXQueue.empty())
	{
		dummySurface = minXQueue.front();
		minXQueue.pop();
		gmshOutput 	<< dummySurface;
		if(!minXQueue.empty())
		{
			// Then the list still isn't empty, add comma
			gmshOutput << ", ";
		}
		else
		{
			// Then the list ended
			gmshOutput << "};" << std::endl;
		}
	}
	gmshOutput << std::endl;

	gmshOutput 	<< "Physical Surface(" << BOUNDARY_ID_MAX_X << ") = {";
	while(!maxXQueue.empty())
	{
		dummySurface = maxXQueue.front();
		maxXQueue.pop();
		gmshOutput 	<< dummySurface;
		if(!maxXQueue.empty())
		{
			// Then the list still isn't empty, add comma
			gmshOutput << ", ";
		}
		else
		{
			// Then the list ended
			gmshOutput << "};" << std::endl;
		}
	}
	gmshOutput << std::endl;

	gmshOutput 	<< "Physical Surface(" << BOUNDARY_ID_MIN_Y << ") = {";
	while(!minYQueue.empty())
	{
		dummySurface = minYQueue.front();
		minYQueue.pop();
		gmshOutput 	<< dummySurface;
		if(!minYQueue.empty())
		{
			// Then the list still isn't empty, add comma
			gmshOutput << ", ";
		}
		else
		{
			// Then the list ended
			gmshOutput << "};" << std::endl;
		}
	}

	gmshOutput 	<< "Physical Surface(" << BOUNDARY_ID_MAX_Y << ") = {";
	while(!maxYQueue.empty())
	{
		dummySurface = maxYQueue.front();
		maxYQueue.pop();
		gmshOutput 	<< dummySurface;
		if(!maxYQueue.empty())
		{
			// Then the list still isn't empty, add comma
			gmshOutput << ", ";
		}
		else
		{
			// Then the list ended
			gmshOutput << "};" << std::endl;
		}
	}
	gmshOutput << std::endl;

	gmshOutput 	<< "Physical Surface(" << BOUNDARY_ID_MIN_Z << ") = {";
	while(!minZQueue.empty())
	{
		dummySurface = minZQueue.front();
		minZQueue.pop();
		gmshOutput 	<< dummySurface;
		if(!minZQueue.empty())
		{
			// Then the list still isn't empty, add comma
			gmshOutput << ", ";
		}
		else
		{
			// Then the list ended
			gmshOutput << "};" << std::endl;
		}
	}
	gmshOutput << std::endl;

	gmshOutput 	<< "Physical Surface(" << BOUNDARY_ID_MAX_Z << ") = {";
	while(!maxZQueue.empty())
	{
		dummySurface = maxZQueue.front();
		maxZQueue.pop();
		gmshOutput 	<< dummySurface;
		if(!maxZQueue.empty())
		{
			// Then the list still isn't empty, add comma
			gmshOutput << ", ";
		}
		else
		{
			// Then the list ended
			gmshOutput << "};" << std::endl;
		}
	}
	gmshOutput << std::endl;

	gmshOutput.close();

	centerOutput.close();
}

void ExportVoronoiToGmshTest(voro::container &con, std::string &outFilename)
{
	std::ofstream gmshOutput(outFilename, std::ofstream::trunc);

	// --> Print the information of each face

	// Vector containing the vertices
	/*
	 * --> The information of a face "k" with "n_k" vertices is located in the
	 *     middle of "face_vertices", with the structure
	 *
	 *     ... n_k, v_{0,n_k}, v_{1,n_k}, ..., v_{n_k - 1,n_k}, ...
	 *
	 *     where v_{i,n_k} is the i'th vertex. The positions of the "n_k"'s
	 *     follow a recursion
	 *
	 *     pos_{n_0} 		= 0
	 *     pos_{n_1} 		= (n_0 + 1)				= pos_{n_0} + (n_0 + 1)
	 *     pos_{n_2} 		= (n_0 + 1) + (n_1 + 1)	= pos_{n_1} + (n_1 + 1)
	 *
	 *     pos_{n_(k+1)}	= 						= pos_{n_k} + (n_k + 1)
	 *
	 *     we will use the auxiliary index "aux_index" to save it
	 */
	std::vector<int>	face_vertices;
	int aux_index;
	int nb_of_vertices;

	// Cell neighbors
	std::vector<int>	cell_neighbours;

	// Cell position
	double pos_x, pos_y, pos_z;

	// Cell ID
	int cell_id;

	// Loop variable used to run over all cells
	voro::c_loop_all	cell_loop(con);

	// Dummy cell with neighbour informations
	voro::voronoicell_neighbor dummy_cell;

	// Dummy boolean used as a marker to when to stop iterating
	bool	keepIterating = cell_loop.start();

	while(keepIterating)
	{
		// Voronoi structure is properly set

		// Transfer needed data to "dummy_cell"
		con.compute_cell(dummy_cell,cell_loop);

		// ... and extract it
		cell_loop.pos(pos_x,pos_y,pos_z);
		cell_id = cell_loop.pid();

		dummy_cell.neighbors(cell_neighbours);
		dummy_cell.face_vertices(face_vertices);

		aux_index = 0;

		gmshOutput 	<< " Cell no. " << cell_id << " @ ("
					<< pos_x << ", " << pos_y << ", " << pos_z << ")"
					<< std::endl;

		for(unsigned int iii = 0; iii < cell_neighbours.size(); ++iii)
		{
			nb_of_vertices = face_vertices[aux_index];
//			if(cell_neighbours[iii] > cell_id || cell_neighbours[iii] < 0)
			{
				// Then we didn't treat this face yet
				gmshOutput	<< " -- Face no. " << iii << ", " << dummy_cell.p << " vertices" << std::endl;
				gmshOutput	<< "             ";
				for(int jjj = 1; jjj < nb_of_vertices + 1; ++jjj)
				{
					gmshOutput << face_vertices[aux_index + jjj] << " ";
				}
				gmshOutput << std::endl;
			}
			// In any case, jump to the next face
			aux_index += nb_of_vertices + 1;
		}

		// Go to next cell (returns "false" if it's over)
		keepIterating = cell_loop.inc();
	}

	gmshOutput.close();
}

int main(int argc, char *argv[])
{
	GetPot cmd_line(argc,argv);

	// Set up constants for the container geometry
	double x_min = -1, x_max = 1;
	double y_min = -1, y_max = 1;
	double z_min = -1, z_max = 1;

	// Set the number of particles that are going to be randomly introduced
	int particles = 20;

	// Set up periodicities
	bool xIsPeriodic = false;
	bool yIsPeriodic = false;
	bool zIsPeriodic = false;

	// Set up weight
	double weight = 1.0;

	// Set up filename
	std::string output_geo_filename("test_voronoi_to_gmsh");

	// --> Read parameters from command line
	// - Filename
	if ( cmd_line.search(1, "-o") )
	{
		output_geo_filename = cmd_line.next(output_geo_filename.c_str());
	}

	// - Dimensions
	if ( cmd_line.search(1, "-xmin") )
	{
		x_min = cmd_line.next(x_min);
	}

	if ( cmd_line.search(1, "-xmax") )
	{
		x_max = cmd_line.next(x_max);
	}

	if ( cmd_line.search(1, "-ymin") )
	{
		y_min = cmd_line.next(y_min);
	}

	if ( cmd_line.search(1, "-ymax") )
	{
		y_max = cmd_line.next(y_max);
	}

	if ( cmd_line.search(1, "-zmin") )
	{
		z_min = cmd_line.next(z_min);
	}

	if ( cmd_line.search(1, "-zmax") )
	{
		z_max = cmd_line.next(z_max);
	}

	// Number of particles!
	if ( cmd_line.search(1, "-p") )
	{
		particles = cmd_line.next(particles);
	}

	// Output filename
	if ( cmd_line.search(1, "-o") )
	{
		output_geo_filename = cmd_line.next(output_geo_filename.c_str());
	}

	if ( cmd_line.search(1, "-w") )
	{
		weight = cmd_line.next(weight);
		if(weight > 1)
		{
			weight = 1.0;
		}
	}

	// Set up the number of blocks that the container is divided into
	int part_per_block = 8;

	double blocks_density = (double) particles / part_per_block;
	double blocks_linear_density = pow(blocks_density,1/3.);
	int nx = 2*blocks_linear_density;
	int ny = 2*blocks_linear_density;
	int nz = 2*blocks_linear_density;

	// Random variables
	boost::random::lagged_fibonacci607 m_rng;
	boost::random::uniform_01<> rnd;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,nx,ny,nz,
			xIsPeriodic,yIsPeriodic,zIsPeriodic,part_per_block);

	double x,y,z;

	// Randomly add particles into the container
	for(int iii=0; iii<particles; iii++) {
		x=x_min+rnd(m_rng)*(x_max-x_min);
		y=y_min+rnd(m_rng)*(y_max-y_min);
		z=z_min+rnd(m_rng)*(z_max-z_min);
		con.put(iii,x,y,z);
	}

	// Export it to Gmsh format
	ExportVoronoiToGmsh(con,weight,output_geo_filename);

	// And build physical parameters
	double youngMean = 200;
	double muMean = 80;

	double youngAmpl = youngMean*0.05;
	double muAmpl = muMean*0.05;

	std::string  outPhysicalParameters = output_geo_filename + "_physical.dat";
	std::ofstream physOutput(outPhysicalParameters, std::ofstream::trunc);

	double youngValue = -1;
	double muValue = -1;

	physOutput << particles << std::endl;
	for(int iii=0; iii<particles; iii++)
	{
		youngValue = youngMean + (2*rnd(m_rng) - 1) * youngAmpl;
		muValue = muMean + (2*rnd(m_rng) - 1) * muAmpl;
		physOutput << youngValue << " " << muValue << " " << iii + 1 << std::endl;
	}
	physOutput.close();
}
