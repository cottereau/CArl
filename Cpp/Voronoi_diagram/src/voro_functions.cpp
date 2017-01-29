/*
 * voro_functions.cpp
 *
 *  Created on: Jan 29, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#include "voro_functions.h"


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

void voro_statistics(voro::container &con, double& nb_of_faces, double& nb_of_edges)
{
	// Loop variable used to run over all cells
	voro::c_loop_all	cell_loop(con);

	// Dummy cell with neighbour informations
	voro::voronoicell dummy_cell;

	// Dummy boolean used as a marker to when to stop iterating
	bool keepIterating = cell_loop.start();
	std::vector<int> dummy_face_edges;

	nb_of_faces = 0;
	nb_of_edges = 0;

	while(keepIterating)
	{
		// Transfer needed data to "dummy_cell"
		con.compute_cell(dummy_cell,cell_loop);

		nb_of_faces += dummy_cell.number_of_faces();

		dummy_cell.face_orders(dummy_face_edges);
		for(unsigned int iii = 0; iii < dummy_face_edges.size(); ++iii)
		{
			nb_of_edges+= dummy_face_edges[iii];
		}

		keepIterating = cell_loop.inc();
	}

	nb_of_edges /= nb_of_faces;
	nb_of_faces /= con.total_particles();
//	std::cout << nb_of_faces / con.total_particles() << " " << nb_of_edges / nb_of_faces << std::endl;
}

void NewExportVoronoiToGmsh(voro::container &con, double weight, std::string &baseFilename)
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

	std::vector<double> point_min = { con.ax, con.ay, con.az };
	std::vector<double> point_max = { con.bx, con.by, con.bz };

	gmshOutput << "// Number of cells : " << con.total_particles() << std::endl
			   << "// Dimensions      : (" << con.ax << ", " << con.bx << ") ("
			                               << con.ay << ", " << con.by << ") ("
										   << con.az << ", " << con.bz << ")"
			   << std::endl
	           << "// Volume          : " << con.sum_cell_volumes()
	           << std::endl << std::endl;

	gmshOutput << "// Weights" << std::endl;
    gmshOutput << "w1 = " << weight << ";" << std::endl << std::endl;

	// -> First step: estimate the number of vertices, edges and faces
	// "Total" numbers
	int nb_of_vertices_preamble = 0;
	int nb_of_edges_preamble = 0;
	int nb_of_faces_preamble = 0;

	// Dummy boolean used as a marker to when to stop iterating
	bool keepIterating = cell_loop.start();
	std::vector<int> dummy_face_edges;

	while(keepIterating)
	{
		// Transfer needed data to "dummy_cell"
		con.compute_cell(dummy_cell,cell_loop);

		nb_of_vertices_preamble += dummy_cell.p;
		nb_of_faces_preamble += dummy_cell.number_of_faces();

		dummy_cell.face_orders(dummy_face_edges);
		for(unsigned int iii = 0; iii < dummy_face_edges.size(); ++iii)
		{
			nb_of_edges_preamble += dummy_face_edges[iii];
		}

		keepIterating = cell_loop.inc();
	}

	// Preallocate the data structures
	std::vector<std::vector<double> > 	points_list(nb_of_vertices_preamble,std::vector<double>(3,0));
	std::vector<std::vector<int> > 		edges_points_list(nb_of_edges_preamble,std::vector<int>(2,0));
	std::vector<std::vector<int> > 		faces_edges_list(nb_of_faces_preamble);
	std::vector<std::vector<int> > 		cells_faces_list(con.total_particles());
	std::vector<int>					nb_of_points_per_edge(nb_of_edges_preamble,2);

	std::unordered_map<std::vector<long>, int, PointHash_3D, PointHash_3D_Equal > discrete_points;
	discrete_points.reserve(nb_of_vertices_preamble);

	std::vector< int > reduced_point_list(nb_of_vertices_preamble);



	// -> Second step: process the Voronoi cells data into a format compatible
	//    with gmsh

	// Vertices !
	std::vector<double> cell_vertices;
	int cell_nb_of_vertices = 0;

	// Edges and faces!
	int treated_vertices = 1;
	int treated_edges    = 1;
	int treated_faces    = 1;

	int cell_number_of_faces = 0;

	int current_vertex	 = 0;

	int edge_vertex_a	 = 0;
	int edge_vertex_b	 = 0;

	int vertex_global_index = 1;
	int edge_global_index = 1;
	int face_global_index = 1;
	int cell_global_index = 1;

	int vertex_full_global_index = 1;

	double pointX = -1;
	double pointY = -1;
	double pointZ = -1;

	double point_aX = -1;
	double point_aY = -1;
	double point_aZ = -1;

	double point_bX = -1;
	double point_bY = -1;
	double point_bZ = -1;

	double edge_length = 0;

	int minXCounter = 0;
	int maxXCounter = 0;

	int minYCounter = 0;
	int maxYCounter = 0;

	int minZCounter = 0;
	int maxZCounter = 0;

	centerOutput << con.total_particles() << std::endl;

	// Dummy boolean used as a marker to when to stop iterating
	keepIterating = cell_loop.start();

	// Methods used in the hash functions
	point_hasher vertex_hash;
	vertex_hash.set_grid_constraints(point_min,point_max);

	std::vector<double> dummy_point(3);
	std::vector<long> dummy_discrete_point(3);

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
			// Hash the point
			for(int jjj = 0; jjj < 3; ++jjj)
			{
				dummy_point[jjj] = cell_vertices[3*iii+jjj];
			}

			vertex_hash.convert_to_discrete(dummy_point,dummy_discrete_point);
			if(discrete_points.find(dummy_discrete_point) == discrete_points.end())
			{
				// New vertex! Add it to the lists
				discrete_points[dummy_discrete_point] = vertex_global_index;
				points_list[vertex_global_index - 1] = dummy_point;
				++vertex_global_index;
			}

			reduced_point_list[vertex_full_global_index] = discrete_points[dummy_discrete_point];
//			for(int jjj = 0; jjj < 3; ++jjj)
//			{
//				points_list[vertex_global_index - 1][jjj] = cell_vertices[3*iii+jjj];
//			}

			++vertex_full_global_index;
		}

		// Get neighbors and face vertices
		dummy_cell.face_vertices(face_vertices_index);

		cell_number_of_faces = dummy_cell.number_of_faces();

		aux_index = 0;

		// Loop over all the faces
		for(int iii = 0; iii < cell_number_of_faces; ++iii)
		{
			face_nb_of_vertices = face_vertices_index[aux_index];

//			// -> Extract the edges
//			for(int jjj = 1; jjj < face_nb_of_vertices; ++jjj)
//			{
//				dummy_edge[0] = treated_vertices + face_vertices_index[aux_index + jjj];
//				dummy_edge[1] = treated_vertices + face_vertices_index[aux_index + jjj + 1];
//				edges_points_list[edge_global_index - 1][0] = treated_vertices + face_vertices_index[aux_index + jjj];
//				edges_points_list[edge_global_index - 1][1] = treated_vertices + face_vertices_index[aux_index + jjj + 1];
//
//				++edge_global_index;
//			}
//
//			edges_points_list[edge_global_index - 1][0] = treated_vertices + face_vertices_index[aux_index + face_nb_of_vertices];
//			edges_points_list[edge_global_index - 1][1] = treated_vertices + face_vertices_index[aux_index + 1];
//
//			++edge_global_index;

			// -> Extract the edges
			for(int jjj = 1; jjj < face_nb_of_vertices; ++jjj)
			{
				edge_vertex_a = face_vertices_index[aux_index + jjj];
				edge_vertex_b = face_vertices_index[aux_index + jjj + 1];

				point_aX = cell_vertices[3*edge_vertex_a];
				point_aY = cell_vertices[3*edge_vertex_a + 1];
				point_aZ = cell_vertices[3*edge_vertex_a + 2];

				point_bX = cell_vertices[3*edge_vertex_b];
				point_bY = cell_vertices[3*edge_vertex_b + 1];
				point_bZ = cell_vertices[3*edge_vertex_b + 2];

				edge_length = std::sqrt(std::pow(point_bX - point_aX,2.) +
										std::pow(point_bY - point_aY,2.) +
										std::pow(point_bZ - point_aZ,2.));

//				if(edge_length > weight)
//				{
					nb_of_points_per_edge[edge_global_index - 1] = std::ceil(edge_length / weight) + 1;
//				}
//				else
//				{
//
//				}

				edges_points_list[edge_global_index - 1][0] = reduced_point_list[treated_vertices + edge_vertex_a ];
				edges_points_list[edge_global_index - 1][1] = reduced_point_list[treated_vertices + edge_vertex_b ];


				++edge_global_index;
			}

			edge_vertex_a = face_vertices_index[aux_index + face_nb_of_vertices];
			edge_vertex_b = face_vertices_index[aux_index + 1];

			point_aX = cell_vertices[3*edge_vertex_a];
			point_aY = cell_vertices[3*edge_vertex_a + 1];
			point_aZ = cell_vertices[3*edge_vertex_a + 2];

			point_bX = cell_vertices[3*edge_vertex_b];
			point_bY = cell_vertices[3*edge_vertex_b + 1];
			point_bZ = cell_vertices[3*edge_vertex_b + 2];

			edge_length = std::sqrt(std::pow(point_bX - point_aX,2.) +
									std::pow(point_bY - point_aY,2.) +
									std::pow(point_bZ - point_aZ,2.));

//				if(edge_length > weight)
//				{
			nb_of_points_per_edge[edge_global_index - 1] = std::ceil(edge_length / weight) + 1;
//				}
//				else
//				{
//
//				}

			edges_points_list[edge_global_index - 1][0] = reduced_point_list[treated_vertices + edge_vertex_a ];
			edges_points_list[edge_global_index - 1][1] = reduced_point_list[treated_vertices + edge_vertex_b ];

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
			for(int jjj = 0; jjj < face_nb_of_vertices - 1; ++jjj)
			{
				faces_edges_list[face_global_index - 1].push_back(treated_edges + jjj);
			}

			// -> Close the face
			faces_edges_list[face_global_index - 1].push_back(treated_edges + face_nb_of_vertices - 1);
			// Update face index
			++face_global_index;
			treated_edges = edge_global_index;

			// -> Jump to the next face
			aux_index += face_nb_of_vertices + 1;
		}

		// -> Finally, build the cell!
		for(int iii = 0; iii < cell_number_of_faces - 1; ++iii)
		{
			cells_faces_list[cell_global_index - 1].push_back(treated_faces + iii);
		}

		// -> Close the cell
		cells_faces_list[cell_global_index - 1].push_back(treated_faces + cell_number_of_faces - 1);

		// -> Print the center
		centerOutput 	<< pos_x << " " << pos_y << " " << pos_z << " "
						<< cell_global_index << std::endl;

		// Update the cell index
		++cell_global_index;
		treated_faces = face_global_index;

		// Go to next cell (returns "false" if it's over)
		treated_vertices += dummy_cell.p;

		// Go to next cell (returns "false" if it's over)
		keepIterating = cell_loop.inc();
	}

	// -> Third step: reduce the points that are too near from each other
	// -> Fourth step: print the geometry data

	// Points!
	std::cout << vertex_global_index << " " << edge_global_index << " " << face_global_index << " " << cell_global_index << std::endl;
	for(int iii = 1; iii < vertex_global_index; ++iii)
	{
		gmshOutput 	<< "Point(" << iii << ") = {";
		for(int jjj = 0; jjj < 2; ++jjj)
		{
			gmshOutput << points_list[iii - 1][jjj] << ", ";
		}
		gmshOutput << points_list[iii - 1][2] << "};" << std::endl;
	}
	gmshOutput << std::endl;

	// Edges!
	for(int iii = 1; iii < edge_global_index; ++iii)
	{
		gmshOutput	<< "Line(" << iii << ") = {"
					<< edges_points_list[iii - 1][0] << ", " << edges_points_list[iii - 1][1] << "};"
					<< std::endl;
		gmshOutput  << "Transfinite Line(" << iii << ") = "
				<< nb_of_points_per_edge[iii - 1] << " Using Bump 0.25;" << std::endl;
	}
	gmshOutput << std::endl;

	// Faces!
	for(int iii = 1; iii < face_global_index; ++iii)
	{
		gmshOutput	<< "Line Loop(" << iii << ") = {";

		for(unsigned int jjj = 0; jjj < faces_edges_list[iii - 1].size() - 1; ++jjj)
		{
			gmshOutput	<< faces_edges_list[iii - 1][jjj] << ", ";
		}

		// -> Close the face
		gmshOutput	<< faces_edges_list[iii - 1][faces_edges_list[iii - 1].size() - 1] << "};" << std::endl;
		gmshOutput	<< "Plane Surface(" <<  iii << ") = {" <<  iii << "};" << std::endl;
//		gmshOutput  << "Transfinite Surface {"  <<  iii << "};" << std::endl;
	}

	// Cells!
	for(int iii = 1; iii < cell_global_index; ++iii)
	{
		gmshOutput << "Surface Loop(" << iii << ") = {";
		for(unsigned int jjj = 0; jjj < cells_faces_list[iii - 1].size() - 1; ++jjj)
		{
			gmshOutput	<< cells_faces_list[iii - 1][jjj] << ", ";
		}

		// -> Close the cell
		gmshOutput	<< cells_faces_list[iii - 1][cells_faces_list[iii - 1].size() - 1] << "};" << std::endl;

		// -> Set up the volumes
		gmshOutput	<< "Volume(" <<  iii
					<< ") = {" <<  iii << "};" << std::endl;
		gmshOutput	<< "Physical Volume(" <<  iii
					<< ") = {" <<  iii << "};" << std::endl;
		gmshOutput	<< std::endl;
	}

	gmshOutput << "Coherence;" << std::endl;

	// -> Fifth step: print the physics data
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
