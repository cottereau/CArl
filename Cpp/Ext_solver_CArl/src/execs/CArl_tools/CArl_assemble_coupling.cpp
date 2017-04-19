#include "CArl_assemble_coupling.h"

libMesh::ExplicitSystem& add_explicit_elasticity(	libMesh::EquationSystems& input_systems,
												libMesh::Order order = libMesh::FIRST,
												libMesh::FEFamily family = libMesh::LAGRANGE)
{
	libMesh::ExplicitSystem& elasticity_system =
			input_systems.add_system<libMesh::ExplicitSystem> ("Elasticity");

	elasticity_system.add_variable("u", order, family);
	elasticity_system.add_variable("v", order, family);
	elasticity_system.add_variable("w", order, family);

	return elasticity_system;
}

libMesh::ImplicitSystem& add_elasticity(	libMesh::EquationSystems& input_systems,
												libMesh::Order order = libMesh::FIRST,
												libMesh::FEFamily family = libMesh::LAGRANGE)
{
	libMesh::ImplicitSystem& elasticity_system =
			input_systems.add_system<libMesh::ImplicitSystem> ("Elasticity");

	elasticity_system.add_variable("u", order, family);
	elasticity_system.add_variable("v", order, family);
	elasticity_system.add_variable("w", order, family);

	return elasticity_system;
}

int main(int argc, char** argv) {

	// --- Initialize libMesh
	libMesh::LibMeshInit init(argc, argv);

	// Do performance log?
	const bool MASTER_bPerfLog_carl_libmesh = true;
	libMesh::PerfLog perf_log("Main program", MASTER_bPerfLog_carl_libmesh);

	// libMesh's C++ / MPI communicator wrapper
	libMesh::Parallel::Communicator& WorldComm = init.comm();

	// Number of processors and processor rank.
	int rank = WorldComm.rank();
	int nodes = WorldComm.size();

	// Create local communicator
	libMesh::Parallel::Communicator LocalComm;
	WorldComm.split(rank,rank,LocalComm);

	// --- Set up inputs

	// Command line parser
	GetPot command_line(argc, argv);

	// File parser
	GetPot field_parser;

	// If there is an input file, parse it to get the parameters. Else, parse the command line
	std::string input_filename;
	if (command_line.search(2, "--inputfile", "-i")) {
		input_filename = command_line.next(input_filename);
		field_parser.parse_input_file(input_filename, "#", "\n", " \t\n");
	} else {
		field_parser = command_line;
	}

	carl::coupling_assemble_coupling_input_params input_params;
	carl::get_assemble_coupling_input_params(field_parser, input_params);

	// Check libMesh installation dimension
	const unsigned int dim = 3;

	libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

	// --- Set up the meshes

	// - Parallelized meshes: A, B, mediator and weight
	perf_log.push("Meshes - Parallel","Read files:");
	libMesh::Mesh mesh_BIG(WorldComm, dim);
	mesh_BIG.read(input_params.mesh_BIG_file);
	mesh_BIG.prepare_for_use();

	libMesh::Mesh mesh_micro(WorldComm, dim);
	mesh_micro.read(input_params.mesh_micro_file);
	mesh_micro.prepare_for_use();

	libMesh::Mesh mesh_mediator(WorldComm, dim);
	mesh_mediator.allow_renumbering(false);
	mesh_mediator.read(input_params.mesh_mediator_file);
	mesh_mediator.prepare_for_use();

	perf_log.pop("Meshes - Parallel","Read files:");

	// - Local meshes: restrict A and restrict B

	// -> libMesh's default mesh IS the ReplicatedMesh, which creates a copy of
	//    itself on each processor, but partitions the iterators over each
	//    processor to allow the parallelization of the code. The usage of
	//    ParallelMesh is still a bit risky, and, performance-wise, is worse
	//    than ReplicatedMesh for less than a few thousand processors.

	perf_log.push("Meshes - Serial","Read files:");
	libMesh::ReplicatedMesh mesh_R_BIG(WorldComm, dim);
	mesh_R_BIG.allow_renumbering(false);
	mesh_R_BIG.read(input_params.mesh_restrict_BIG_file);
	mesh_R_BIG.prepare_for_use();

	libMesh::ReplicatedMesh mesh_R_micro(WorldComm, dim);
	mesh_R_micro.allow_renumbering(false);
	mesh_R_micro.read(input_params.mesh_restrict_micro_file);
	mesh_R_micro.prepare_for_use();

	// - Local mesh: intersection mesh
	libMesh::Mesh mesh_inter(LocalComm, dim);
	mesh_inter.allow_renumbering(false);
	std::string local_inter_mesh_filename = input_params.common_inter_file + "_r_"
										+ std::to_string(rank) + "_n_" + std::to_string(nodes) + ".e";
	std::string local_inter_table_filename = input_params.common_inter_file + "_r_"
										+ std::to_string(rank) + "_n_" + std::to_string(nodes) + "_inter_table_Full.dat";
	std::string global_inter_table_filename = input_params.common_inter_file + "_stitched_inter_table_Full.dat";

	mesh_inter.read(local_inter_mesh_filename);
	mesh_inter.prepare_for_use();

	perf_log.pop("Meshes - Serial","Read files:");

	// --- Prepare the equivalence tables and the intersection mappings
	perf_log.push("Equivalence / intersection tables","Read files:");

	// Local intersection pairs (system meshes)
	std::unordered_map<int,std::pair<int,int> > local_intersection_pairs_map;

	// Local intersection pairs (restricted meshes)
	std::unordered_map<int,std::pair<int,int> > local_intersection_restricted_pairs_map;

	// Intersection elements to be treated by this processor
	std::unordered_map<int,int> local_intersection_meshI_to_inter_map;

	// Equivalence tables between system and restricted meshes (and vice-versa)
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

	// Set the local intersection tables
	if(input_params.mediator_type == carl::MediatorType::USE_MACRO)
	{
		carl::set_local_intersection_tables(
				WorldComm,
				mesh_inter,
				local_inter_table_filename,
				input_params.equivalence_table_restrict_BIG_file,
				input_params.equivalence_table_restrict_micro_file,

				equivalence_table_BIG_to_R_BIG,
				equivalence_table_micro_to_R_micro,

				local_intersection_pairs_map,
				local_intersection_restricted_pairs_map,
				local_intersection_meshI_to_inter_map);
	}
	else if(input_params.mediator_type == carl::MediatorType::USE_MICRO)
	{
		carl::set_local_intersection_tables(
				WorldComm,
				mesh_inter,
				local_inter_table_filename,
				input_params.equivalence_table_restrict_micro_file,
				input_params.equivalence_table_restrict_BIG_file,

				equivalence_table_micro_to_R_micro,
				equivalence_table_BIG_to_R_BIG,

				local_intersection_pairs_map,
				local_intersection_restricted_pairs_map,
				local_intersection_meshI_to_inter_map);
	}

	// Build mappings with the intersections of the mediator mesh with the micro and macro meshes - used to preallocate the coupling matrices
	std::unordered_multimap<int,int> inter_mediator_BIG;
	std::unordered_multimap<int,int> inter_mediator_micro;
	carl::set_global_mediator_system_intersection_lists(
			WorldComm,
			global_inter_table_filename,
			equivalence_table_BIG_to_R_BIG,
			equivalence_table_R_BIG_to_BIG,
			inter_mediator_BIG,
			inter_mediator_micro);
	perf_log.pop("Equivalence / intersection tables","Read files:");

	// --- Generate the equation systems
	perf_log.push("Initialization","System initialization:");
	carl::assemble_coupling_matrices CoupledMatrices(WorldComm);

	libMesh::EquationSystems& equation_systems_inter =
					CoupledMatrices.add_inter_EquationSystem("InterSys", mesh_inter);

	libMesh::ExplicitSystem& elasticity_system_inter
										= add_explicit_elasticity(equation_systems_inter);

	carl::reduced_system_init(elasticity_system_inter);

	perf_log.pop("Initialization","System initialization:");

	// - Build the MACRO system

	perf_log.push("Macro system","System initialization:");
	libMesh::EquationSystems& equation_systems_BIG =
					CoupledMatrices.set_BIG_EquationSystem("BigSys", mesh_BIG);

	// [MACRO] Set up the physical properties
	libMesh::ExplicitSystem& elasticity_system_BIG
										= add_explicit_elasticity(equation_systems_BIG);

	carl::reduced_system_init(elasticity_system_BIG);
	perf_log.pop("Macro system","System initialization:");

	// - Build the micro system

	perf_log.push("Micro system","System initialization:");
	libMesh::EquationSystems& equation_systems_micro =
					CoupledMatrices.add_micro_EquationSystem<libMesh::PetscMatrix<libMesh::Number> >("MicroSys", mesh_micro);

	// [MICRO] Set up the physical properties
	libMesh::ExplicitSystem& elasticity_system_micro
										= add_explicit_elasticity(equation_systems_micro);

	carl::reduced_system_init(elasticity_system_micro);
	perf_log.pop("Micro system","System initialization:");

	// - Build the RESTRICTED BIG system
	perf_log.push("RESTRICTED macro system","System initialization:");
	libMesh::EquationSystems& equation_systems_R_BIG =
					CoupledMatrices.set_Restricted_BIG_EquationSystem("BigSys", mesh_R_BIG);

	// [R. MACRO] Set up the physical properties
	libMesh::ExplicitSystem& elasticity_system_R_BIG
										= add_explicit_elasticity(equation_systems_R_BIG);

	carl::reduced_system_init(elasticity_system_R_BIG);

	perf_log.pop("RESTRICTED macro system","System initialization:");

	// - Build the RESTRICTED micro system

	perf_log.push("RESTRICTED micro system","System initialization:");
	libMesh::EquationSystems& equation_systems_R_micro =
					CoupledMatrices.add_Restricted_micro_EquationSystem("MicroSys", mesh_R_micro);

	// [R. MICRO] Set up the physical properties
	libMesh::ExplicitSystem& elasticity_system_R_micro
										= add_explicit_elasticity(equation_systems_R_micro);

	carl::reduced_system_init(elasticity_system_R_micro);
	perf_log.pop("RESTRICTED micro system","System initialization:");

	// - Build the mediator system

	perf_log.push("Mediator system","System initialization:");
	libMesh::EquationSystems& equation_systems_mediator =
					CoupledMatrices.add_mediator_EquationSystem("MediatorSys", mesh_mediator);

	libMesh::ExplicitSystem& elasticity_system_mediator
										= add_explicit_elasticity(equation_systems_mediator);

	carl::reduced_system_init(elasticity_system_mediator);

	perf_log.pop("Mediator system","System initialization:");

	// Set parameters
	CoupledMatrices.set_coupling_parameters("MicroSys",input_params.coupling_rigidity,input_params.coupling_width);

	// Use H1 norm
	CoupledMatrices.use_H1_coupling("MicroSys");

	// Assemble!
	CoupledMatrices.assemble_coupling_elasticity_3D_parallel("BigSys","MicroSys",
			"InterSys","MediatorSys",
			mesh_R_BIG, mesh_R_micro,
			local_intersection_pairs_map,
			local_intersection_restricted_pairs_map,
			local_intersection_meshI_to_inter_map,
			inter_mediator_BIG,
			inter_mediator_micro);

	// Print matrices!
	CoupledMatrices.print_matrices_matlab("MicroSys",input_params.output_base);
	CoupledMatrices.print_PETSC_matrices("MicroSys",input_params.output_base);

	return 0;
}
