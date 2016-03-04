#include "main.h"

/*
 * 		Sketch of the first version of the program -----------------------------
 *
 * 		-> Read the meshes A, B and the coupling region file.
 *
 * 		-> Read the intersection mesh I.
 *
 * 		-> Read the mediator mesh, M. (!)
 *
 * 		-> Read the mesh restrictions, R_A and R_B, and the equivalence tables
 * 		   between them and the original meshes, t_R_A->A and t_R_B->B.
 *
 * 		-> Copy R_A, R_B and I to all processors. (!)
 *
 * 		-> Run the coupling assemble program using the mediator space as the
 * 		   main loop index.
 *
 * 		-> Use the equivalence tables to fill the coupling matrices.
 *
 * 		-> Export them using PETSc MatView methods.
 *
 * 		Optimizations ----------------------------------------------------------
 *
 * 		-> Optimize the partitioning of M using MeshTools: build a weighting
 * 		   vector for the intersections, w_I, convert it to a w_M weighting, and
 * 		   partition M.
 *
 * 		-> Partition I by hand using this information.
 *
 * 		-> Partition R_A and R_B using this information.
 *
 * 		Functions I'll need ----------------------------------------------------
 *
 * 		TODO :	read meshes A, B and I in a way to make a single copy on each
 * 				processor - watch out for collisions!
 *
 * 		TODO :	read the intersection and the restriction equivalence tables.
 * 				Divide them on processor-by-processor tables.
 *
 * 		TODO :	adapt the coupling assemble methods in a way such that the main
 * 				loop is ran over the mediator mesh, and not the intersection
 * 				mesh.
 *
 * 		TODO :	export the matrix using Petsc interfaces.
 *
 * 		TODO :	build a matrix import function.
 *
 * 		Things to watch out for ------------------------------------------------
 *
 * 		-> Matrix renumerotations during the prepare_for_use step:
 * 		   This is done before the partitioning, and, according to the
 * 		   documentation, is done to guarantee that the elements are organized
 * 		   in contiguous blocks. From this, and from the current tests, it
 * 		   SHOULD not be a problem ... but still, watch out ...
 */

struct carl_coupling_generation_input_params {
	std::string physical_params_file;

	std::string mesh_BIG_file;
	std::string mesh_micro_file;

	std::string mesh_restrict_BIG_file;
	std::string mesh_restrict_micro_file;

	std::string mesh_mediator_file;
	std::string mesh_inter_file;

	std::string equivalence_table_restrict_BIG_file;
	std::string equivalence_table_restrict_micro_file;
	std::string equivalence_table_mediator;

	std::string intersection_table_BIG_micro_file;
	std::string intersection_table_I_file;
	std::string intersection_table_full;

	bool b_UseMesh_BIG_AsMediator;
	bool b_UseMesh_micro_AsMediator;
	bool b_UseMesh_extra_AsMediator;

	double meanE;
	double meanMu;
	double coupling_const;
	double mean_distance;

	std::string output_coupling_matrix_med_BIG;
	std::string output_coupling_matrix_med_micro;
};

void get_input_params(GetPot& field_parser,
		carl_coupling_generation_input_params& input_params) {

	// Set mesh files
	if (field_parser.search(3, "--meshA", "-mA", "MeshA")) {
		input_params.mesh_BIG_file = field_parser.next(
				input_params.mesh_BIG_file);
	} else {
		homemade_error_msg("Missing the A mesh file!");
	}

	if (field_parser.search(3, "--meshB", "-mB", "MeshB")) {
		input_params.mesh_micro_file = field_parser.next(
				input_params.mesh_micro_file);
	} else {
		homemade_error_msg("Missing the B mesh file!");
	}

	if (field_parser.search(3, "--meshAR", "-mAR", "Mesh_A_Restriction")) {
		input_params.mesh_restrict_BIG_file = field_parser.next(
				input_params.mesh_restrict_BIG_file);
	} else {
		homemade_error_msg("Missing the restricted A mesh file!");
	}

	if (field_parser.search(3, "--meshBR", "-mBR", "Mesh_B_Restriction")) {
		input_params.mesh_restrict_micro_file = field_parser.next(
				input_params.mesh_restrict_micro_file);
	} else {
		homemade_error_msg("Missing the restricted B mesh file!");
	}

	if (field_parser.search(3, "--meshI", "-mI", "MeshInter")) {
		input_params.mesh_inter_file = field_parser.next(
				input_params.mesh_inter_file);
	} else {
		homemade_error_msg("Missing the intersection mesh file!");
	}

	// Set the equivalence and intersection tables
	if (field_parser.search(2, "--tableRA", "Mesh_A_RestrictionEquivalenceTable")) {
		input_params.equivalence_table_restrict_BIG_file = field_parser.next(
				input_params.equivalence_table_restrict_BIG_file);
	} else {
		homemade_error_msg("Missing the equivalence table for the mesh A!");
	}

	if (field_parser.search(2, "--tableRB", "Mesh_B_RestrictionEquivalenceTable")) {
		input_params.equivalence_table_restrict_micro_file = field_parser.next(
				input_params.equivalence_table_restrict_micro_file);
	} else {
		homemade_error_msg("Missing the equivalence table for the mesh B!");
	}

	if (field_parser.search(2, "--tablePairs", "IntersectionPairsTable")) {
		input_params.intersection_table_BIG_micro_file = field_parser.next(
				input_params.intersection_table_BIG_micro_file);
	} else {
		homemade_error_msg("Missing the intersection pairs file!");
	}

	if (field_parser.search(2, "--tableFullI", "FullIntersectionElementsTable")) {
		input_params.intersection_table_full = field_parser.next(
				input_params.intersection_table_full);
	} else {
		homemade_error_msg("Missing the full intersection elements file!");
	}

	// Set the mediator mesh
	input_params.b_UseMesh_BIG_AsMediator = false;
	input_params.b_UseMesh_micro_AsMediator = false;
	input_params.b_UseMesh_extra_AsMediator = false;
	if (field_parser.search(1,"Use_A_AsMediator"))
	{
		input_params.b_UseMesh_BIG_AsMediator = true;
	}
	if (field_parser.search(1,"Use_B_AsMediator"))
	{
		input_params.b_UseMesh_micro_AsMediator = true;
	}
	if (field_parser.search(1,"Use_extra_AsMediator"))
	{
		input_params.b_UseMesh_extra_AsMediator = true;
	}
	if(input_params.b_UseMesh_BIG_AsMediator
			+ input_params.b_UseMesh_micro_AsMediator
			+ input_params.b_UseMesh_extra_AsMediator > 1)
	{
		homemade_error_msg("Choose only one mesh as mediator!");
	}

	if(input_params.b_UseMesh_BIG_AsMediator)
	{
		input_params.mesh_mediator_file = input_params.mesh_restrict_BIG_file;
		input_params.equivalence_table_mediator = input_params.equivalence_table_restrict_BIG_file;
	}
	if(input_params.b_UseMesh_micro_AsMediator)
	{
		input_params.mesh_mediator_file = input_params.mesh_restrict_micro_file;
		input_params.equivalence_table_mediator = input_params.equivalence_table_restrict_micro_file;
	}
	if(input_params.b_UseMesh_extra_AsMediator)
	{
//		if (field_parser.search(3, "--meshM", "-mM", "MeshMediator")) {
//			input_params.mesh_mediator_file = field_parser.next(
//					input_params.mesh_mediator_file);
//		} else {
//			homemade_error_msg("Missing the mediator mesh file!");
//		}
		libmesh_not_implemented_msg("Still implementing the external mesh case!");
	}

	// Set constant parameters
	if (field_parser.search(3, "-p", "--parameters", "PhysicalParameters")) {
		input_params.physical_params_file = field_parser.next(
				input_params.physical_params_file);
		calculate_average_params(input_params.physical_params_file,
				input_params.meanE,input_params.meanMu);
		input_params.coupling_const = eval_lambda_1(input_params.meanE,
				input_params.meanMu);
	}
	else if (field_parser.search(2, "--meanE", "MeanE") &&
			field_parser.search(2, "--meanMu", "MeanMu")) {
		field_parser.search(2, "--meanE", "MeanE");
		input_params.meanE = field_parser.next(input_params.meanE);
		field_parser.search(2, "--meanMu", "MeanMu");
		input_params.meanMu = field_parser.next(input_params.meanMu);
		input_params.coupling_const = eval_lambda_1(input_params.meanE,
				input_params.meanMu);
	}
	else if (field_parser.search(2, "--coupling", "CouplingConstant"))
	{
		input_params.coupling_const = field_parser.next(
				input_params.coupling_const);
	}
	else
	{
		homemade_error_msg("Missing the coupling constant");
	}

	input_params.mean_distance = 0.2;

	if (field_parser.search(2, "--dist", "CouplingMeshScale")) {
		input_params.mean_distance = field_parser.next(
				input_params.mean_distance);
	}

	// Set output files
	input_params.output_coupling_matrix_med_BIG =
			"meshes/parallel_test/output/coupling_matrix_mediator_A.dat";
	input_params.output_coupling_matrix_med_micro =
			"meshes/parallel_test/output/coupling_matrix_mediator_B.dat";
	if (field_parser.search(3, "-oA", "--outputA", "OutputMatrixA")) {
		input_params.output_coupling_matrix_med_BIG = field_parser.next(
				input_params.output_coupling_matrix_med_BIG);
	}
	if (field_parser.search(3, "-oB", "--outputB", "OutputMatrixB")) {
		input_params.output_coupling_matrix_med_micro = field_parser.next(
				input_params.output_coupling_matrix_med_micro);
	}

}
;

int main(int argc, char** argv) {

	// - Start libmesh ---------------------------------------------------------
	const bool MASTER_bPerfLog_carl_libmesh = true;
	libMesh::LibMeshInit init(argc, argv);

	libMesh::PerfLog perf_log("Main program", MASTER_bPerfLog_carl_libmesh);

	// - Set up inputs
	GetPot command_line(argc, argv);
	GetPot field_parser;
	std::string input_filename;

	if (command_line.search(1, "--inputfile")) {
		input_filename = command_line.next(input_filename);
		field_parser.parse_input_file(input_filename, "#", "\n", " \t\n");
	} else {
		field_parser = command_line;
	}

	carl_coupling_generation_input_params input_params;
	get_input_params(field_parser, input_params);

	const unsigned int dim = 3;

	libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

	// Set up the communicator and the rank variables
	libMesh::Parallel::Communicator& WorldComm = init.comm();
	libMesh::Parallel::Communicator LocalComm;

	int rank = WorldComm.rank();
	int nodes = WorldComm.size();

	WorldComm.split(rank,rank,LocalComm);

	// - Read the meshes -------------------------------------------------------

	// - Global meshes: A, B and mediator
	libMesh::Mesh mesh_BIG(WorldComm, dim);
	std::unordered_map<int,int> mesh_BIG_NodeMap;
	std::unordered_map<int,int> mesh_BIG_ElemMap;
	carl::set_mesh_Gmsh(mesh_BIG,input_params.mesh_BIG_file);
	carl::create_mesh_map(input_params.mesh_BIG_file,
			mesh_BIG_NodeMap,mesh_BIG_ElemMap,WorldComm);

	libMesh::Mesh mesh_micro(WorldComm, dim);
	std::unordered_map<int,int> mesh_micro_NodeMap;
	std::unordered_map<int,int> mesh_micro_ElemMap;
	carl::set_mesh_Gmsh(mesh_micro,input_params.mesh_micro_file);
	carl::create_mesh_map(input_params.mesh_micro_file,
			mesh_micro_NodeMap,mesh_micro_ElemMap,WorldComm);

	libMesh::Mesh mesh_inter(WorldComm, dim);
	std::unordered_map<int,int> mesh_inter_NodeMap_table_to_libmesh;
	std::unordered_map<int,int> mesh_inter_NodeMap_libmesh_to_table;
	std::unordered_map<int,int> mesh_inter_ElemMap_table_to_libmesh;
	std::unordered_map<int,int> mesh_inter_ElemMap_libmesh_to_table;

	carl::set_mesh_Gmsh(mesh_inter,input_params.mesh_inter_file);
	carl::create_mesh_map(input_params.mesh_inter_file,
			mesh_inter_NodeMap_table_to_libmesh,mesh_inter_NodeMap_libmesh_to_table,
			mesh_inter_ElemMap_table_to_libmesh,mesh_inter_ElemMap_libmesh_to_table,WorldComm);

	// - Test: print info per proc
	std::ofstream mesh_info_ofstream;
	mesh_info_ofstream.open("meshes/parallel_test/output/mesh_A_" + std::to_string(rank) + "_info.txt");
	mesh_BIG.print_info(mesh_info_ofstream);
	mesh_info_ofstream.close();

	WorldComm.barrier();

	mesh_info_ofstream.open("meshes/parallel_test/output/mesh_B_" + std::to_string(rank) + "_info.txt");
	mesh_micro.print_info(mesh_info_ofstream);
	mesh_info_ofstream.close();

	WorldComm.barrier();

	mesh_info_ofstream.open("meshes/parallel_test/output/mesh_inter_" + std::to_string(rank) + "_info.txt");
	mesh_inter.print_info(mesh_info_ofstream);
	mesh_info_ofstream.close();

	WorldComm.barrier();

	// - Local meshes: restrict A, restrict B and (gasp!) intersection

	// -> TODO : 	create decent mesh reader. For now, each processor reads the
	//				mesh independently ... (aïe ...)

	libMesh::Mesh mesh_R_BIG(LocalComm, dim);
	std::unordered_map<int,int> mesh_R_BIG_NodeMap;
	std::unordered_map<int,int> mesh_R_BIG_ElemMap;
	carl::set_mesh_Gmsh(mesh_R_BIG,input_params.mesh_restrict_BIG_file);
	carl::create_mesh_map(input_params.mesh_restrict_BIG_file,
			mesh_R_BIG_NodeMap,mesh_R_BIG_ElemMap,WorldComm);

	libMesh::Mesh mesh_R_micro(LocalComm, dim);
	std::unordered_map<int,int> mesh_R_micro_NodeMap;
	std::unordered_map<int,int> mesh_R_micro_ElemMap;
	carl::set_mesh_Gmsh(mesh_R_micro,input_params.mesh_restrict_micro_file);
	carl::create_mesh_map(input_params.mesh_restrict_micro_file,
			mesh_R_micro_NodeMap,mesh_R_micro_ElemMap,WorldComm);

	libMesh::Mesh mesh_mediator(LocalComm, dim);
	std::unordered_map<int,int> mesh_mediator_NodeMap;
	std::unordered_map<int,int> mesh_mediator_ElemMap;
	carl::set_mesh_Gmsh(mesh_mediator,input_params.mesh_mediator_file);
	carl::create_mesh_map(input_params.mesh_mediator_file,
			mesh_mediator_NodeMap,mesh_mediator_ElemMap,WorldComm);

	// - Test: print info per proc
	mesh_info_ofstream.open("meshes/parallel_test/output/mesh_RA_" + std::to_string(rank) + "_info.txt");
	mesh_R_BIG.print_info(mesh_info_ofstream);
	mesh_info_ofstream.close();

	WorldComm.barrier();

	mesh_info_ofstream.open("meshes/parallel_test/output/mesh_RB_" + std::to_string(rank) + "_info.txt");
	mesh_R_micro.print_info(mesh_info_ofstream);
	mesh_info_ofstream.close();

	WorldComm.barrier();

	mesh_info_ofstream.open("meshes/parallel_test/output/mesh_mediator_" + std::to_string(rank) + "_info.txt");
	mesh_mediator.print_info(mesh_info_ofstream);
	mesh_info_ofstream.close();

	std::ofstream mesh_data;
	mesh_data.open("meshes/parallel_test/output/mesh_A_data_" + std::to_string(rank)  + ".dat");
	libMesh::MeshBase::const_element_iterator       el     = mesh_BIG.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el = mesh_BIG.active_local_elements_end();

	for ( ; el != end_el; ++el)
	{
		const libMesh::Elem* elem = *el;
		mesh_data << elem->id() << " " << elem->point(0) << " " << elem->point(1) << " " << elem->point(2)<< " " << elem->point(3)<< std::endl;
	}
	mesh_data.close();

	/* 		To do on the first proc
	 * 		- Read and build the equivalence tables between R_X and X, e_X - DONE
	 * 		- Read and build the intersection pairs table between A and B, p_AB - DONE
	 *		- Read and build the intersection indexes table, I_F - DONE
	 *
	 * 		Broadcast e_X, p_AB, I_F - DONE
	 *
	 * 		Build on each proc
	 * 		- The restricted intersection pairs table, p_R,AB- DONE
	 * 		- A local intersection indexes table, I_L - DONE
	 *
	 * 		Convert the pairs table to the libMesh indexing - TODO
	 */

	std::unordered_map<int,std::pair<int,int> > full_intersection_pairs_map;
	std::unordered_map<int,std::pair<int,int> > full_intersection_restricted_pairs_map;
	std::unordered_map<int,int> local_intersection_meshI_to_inter_map;

	std::unordered_map<int,int> equivalence_table_BIG_to_R_BIG;
	std::unordered_map<int,int> equivalence_table_micro_to_R_micro;
	std::unordered_map<int,int> equivalence_table_R_BIG_to_BIG;
	std::unordered_map<int,int> equivalence_table_R_micro_to_micro;

	//	Start by reading and broadcasting the equivalence tables
	carl::set_equivalence_tables(
			WorldComm,
			input_params.equivalence_table_restrict_BIG_file,
			input_params.equivalence_table_restrict_micro_file,
			equivalence_table_BIG_to_R_BIG,
			equivalence_table_micro_to_R_micro,
			equivalence_table_R_BIG_to_BIG,
			equivalence_table_R_micro_to_micro);

	if(input_params.b_UseMesh_BIG_AsMediator)
	{
		carl::set_intersection_tables(
				WorldComm,
				mesh_inter,
				input_params.intersection_table_full,
				input_params.equivalence_table_restrict_BIG_file,
				input_params.equivalence_table_restrict_micro_file,
				mesh_inter_ElemMap_libmesh_to_table,

				equivalence_table_BIG_to_R_BIG,
				equivalence_table_micro_to_R_micro,

				full_intersection_pairs_map,
				full_intersection_restricted_pairs_map,
				local_intersection_meshI_to_inter_map);
	}
	else if(input_params.b_UseMesh_micro_AsMediator)
	{
		carl::set_intersection_tables(
				WorldComm,
				mesh_inter,
				input_params.intersection_table_full,
				input_params.equivalence_table_restrict_micro_file,
				input_params.equivalence_table_restrict_BIG_file,
				mesh_inter_ElemMap_libmesh_to_table,

				equivalence_table_micro_to_R_micro,
				equivalence_table_BIG_to_R_BIG,

				full_intersection_pairs_map,
				full_intersection_restricted_pairs_map,
				local_intersection_meshI_to_inter_map);
	}

	// TODO Convert the tables to libMesh's notation
	// TODO Parallelize it!
	// For now only doing this for the restricted stuff
//	carl::convert_equivalence_to_libmesh(mesh_BIG_ElemMap,mesh_R_BIG_ElemMap,
//			equivalence_table_BIG_to_R_BIG,equivalence_table_R_BIG_to_BIG);
//
//	carl::convert_equivalence_to_libmesh(mesh_micro_ElemMap,mesh_R_micro_ElemMap,
//			equivalence_table_micro_to_R_micro,equivalence_table_R_micro_to_micro);
//
//	if(input_params.b_UseMesh_BIG_AsMediator)
//	{
//		carl::convert_pairs_to_libmesh(mesh_R_BIG_ElemMap,mesh_R_micro_ElemMap
//				full_intersection_restricted_pairs_map);
//	}
//	else if(input_params.b_UseMesh_micro_AsMediator)
//	{
//		carl::convert_pairs_to_libmesh(mesh_R_micro_ElemMap,mesh_R_BIG_ElemMap
//				full_intersection_restricted_pairs_map);
//	}

	return 0;
}
