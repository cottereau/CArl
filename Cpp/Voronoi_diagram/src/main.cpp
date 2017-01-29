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
#include "voro_functions.h"

int main(int argc, char *argv[])
{
	GetPot cmd_line(argc,argv);

	// Set up constants for the container geometry
	double x_min = -1, x_max = 1;
	double y_min = -1, y_max = 1;
	double z_min = -1, z_max = 1;

	// Set the minimal point distance
	double min_distance = 0.1;

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

	if ( cmd_line.search(1, "-dist") )
	{
		min_distance = cmd_line.next(min_distance);
	}
	else
	{
		min_distance = 0.5*weight;
	}

	bool bSkipMesh = false;
	if ( cmd_line.search(1, "-skipmesh") )
	{
		bSkipMesh = true;
	}

	// Random variables
	boost::random::lagged_fibonacci607 m_rng;
	boost::random::uniform_01<> rnd;

	// Set up the number of blocks that the container is divided into
	int part_per_block = 8;

	double blocks_density = (double) particles / part_per_block;
	double blocks_linear_density = pow(blocks_density,1/3.);
	int nx = 2*blocks_linear_density;
	int ny = 2*blocks_linear_density;
	int nz = 2*blocks_linear_density;

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

	if(!bSkipMesh)
	{
		// Export it to Gmsh format
		NewExportVoronoiToGmsh(con,weight,output_geo_filename);
	}

	// Build heterogeneous physical parameters
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

	// Build an anisotropy file too
	std::string  outAnglesParameters = output_geo_filename + "_angles.dat";
	std::ofstream anglesOutput(outAnglesParameters, std::ofstream::trunc);

	double c11 = 198;
	double c12 = 125;
	double c44 = 122;
	anglesOutput 	<< particles << " "
					<< c11 << " "
					<< c12 << " "
					<< c44 << " "
					<< youngMean << " "
					<< muMean << std::endl;

	double angle_x = 0;
	double angle_y = 0;
	double angle_z = 0;

	for(int iii=0; iii<particles; iii++)
	{
		angle_x = 2*rnd(m_rng)*M_PI;
		angle_y =   rnd(m_rng)*M_PI;
		angle_z = 2*rnd(m_rng)*M_PI;
		anglesOutput << angle_x << " " << angle_y << " " << angle_z << " " << iii + 1 << std::endl;
	}
	anglesOutput.close();
}
