#include "mesh_tables.h"

// Mesh and weight parameter tables generation
void carl::set_weight_function_domain_idx(	std::string &filename,
										int& domain_Idx_BIG,
										int& nb_of_domain_Idx,
										std::vector<int>& domain_Idx_micro,
										std::vector<int>& domain_Idx_coupling
										)
{
	std::ifstream dataF(filename);

	// Buffer string
	std::string bufferLine;
	std::stringstream	dataBuffer;

	// Read info until the file ends
	while(std::getline(dataF,bufferLine))
	{
		if(bufferLine.compare("$MacroDomainIdx")==0)
		{
			dataF >> domain_Idx_BIG;
		}

		if(bufferLine.compare("$NbOfMicroDomainIdx")==0)
		{
			dataF >> nb_of_domain_Idx;

			domain_Idx_micro.resize(nb_of_domain_Idx);
			domain_Idx_coupling.resize(nb_of_domain_Idx);
		}

		if(bufferLine.compare("$MicroDomainIdxs")==0)
		{
			for(int iii = 0; iii < nb_of_domain_Idx; ++iii )
			{
				dataF >> domain_Idx_micro[iii];
			}
		}

		if(bufferLine.compare("$CouplingDomainIdxs")==0)
		{
			for(int iii = 0; iii < nb_of_domain_Idx; ++iii )
			{
				dataF >> domain_Idx_coupling[iii];
			}
		}
	}

	dataF.close();
}

//void carl::create_mesh_map(
//		const std::string &filename,
//		std::unordered_map<int,int> &node_map,
//		std::unordered_map<int,int> &element_map)
//{
//	/*
//	 *  libmesh : nodes and elements start at 0, and are continuous
//	 *  GMSH    : nodes and elements can be discontinuous, libmesh builds a map
//	 *  Abaqus  : not sure, but seems similar to GMSH. libmesh builds a map, but
//	 *            differently from the method used for GMSH
//	 */
//
//	//if(boost::filesystem::path(filename).extension().string().compare(".msh")==0)
//	//{
//		build_mesh_map_Gmsh(filename,node_map,element_map);
//	//}
//	//else
//	//{
//	//	libmesh_error_msg("Error: unknown filetype");
//	//}
//}
//
//void carl::create_mesh_map(
//		const std::string &filename,
//		std::unordered_map<int,int> &node_gmsh_to_libmesh_map,
//		std::unordered_map<int,int> &node_libmesh_to_gmsh_map,
//		std::unordered_map<int,int> &element_gmsh_to_libmesh_map,
//		std::unordered_map<int,int> &element_libmesh_to_gmsh_map
//		)
//{
//	/*
//	 *  libmesh : nodes and elements start at 0, and are continuous
//	 *  GMSH    : nodes and elements can be discontinuous, libmesh builds a map
//	 *  Abaqus  : not sure, but seems similar to GMSH. libmesh builds a map, but
//	 *            differently from the method used for GMSH
//	 */
//
//	//if(boost::filesystem::path(filename).extension().string().compare(".msh")==0)
//	//{
//		build_mesh_map_Gmsh(filename,node_gmsh_to_libmesh_map,node_libmesh_to_gmsh_map,
//				element_gmsh_to_libmesh_map,element_libmesh_to_gmsh_map);
//	//}
//	//else
//	//{
//	//	libmesh_error_msg("Error: unknown filetype");
//	//}
//};
//
//void carl::create_mesh_map(
//		const std::string &filename,
//		std::unordered_map<int,int> &node_map,
//		std::unordered_map<int,int> &element_map,
//		const libMesh::Parallel::Communicator& MeshComm)
//{
//	/*
//	 *  libmesh : nodes and elements start at 0, and are continuous
//	 *  Abaqus  : not sure, but seems similar to GMSH. libmesh builds a map, but
//	 *            differently from the method used for GMSH
//	 */
//
//	int rank = MeshComm.rank();
//
//	if(rank == 0)
//	{
//		if( filename.rfind(".msh") < filename.size() )
//		{
//			// Gmsh file : nodes and elements can be discontinuous, libmesh
//			//             builds a map, we must do the same
//			build_mesh_map_Gmsh(filename,node_map,element_map);
//		}
//		else if( filename.rfind(".e") < filename.size() )
//		{
//			// Covers both ".e" and ".exd" formats
//			// Exodus_II file :
//		}
//		else
//		{
//			homemade_error_msg("Error: unknown filetype");
//		}
//	}
//
//	carl::broadcast_index_unordered_map(node_map,MeshComm);
//	carl::broadcast_index_unordered_map(element_map,MeshComm);
//}
//
//void carl::create_mesh_map(
//		const std::string &filename,
//		std::unordered_map<int,int> &node_gmsh_to_libmesh_map,
//		std::unordered_map<int,int> &node_libmesh_to_gmsh_map,
//		std::unordered_map<int,int> &element_gmsh_to_libmesh_map,
//		std::unordered_map<int,int> &element_libmesh_to_gmsh_map,
//		const libMesh::Parallel::Communicator& MeshComm
//		)
//{
//	/*
//	 *  libmesh : nodes and elements start at 0, and are continuous
//	 *  GMSH    : nodes and elements can be discontinuous, libmesh builds a map
//	 *  Abaqus  : not sure, but seems similar to GMSH. libmesh builds a map, but
//	 *            differently from the method used for GMSH
//	 */
//
//	//if(boost::filesystem::path(filename).extension().string().compare(".msh")==0)
//	//{
//		int rank = MeshComm.rank();
//
//		if(rank == 0)
//		{
//			build_mesh_map_Gmsh(filename,node_gmsh_to_libmesh_map,element_gmsh_to_libmesh_map);
//		}
//		MeshComm.barrier();
//		carl::broadcast_index_unordered_map(node_gmsh_to_libmesh_map,MeshComm);
//		carl::broadcast_index_unordered_map(element_gmsh_to_libmesh_map,MeshComm);
//
//		carl::invert_index_unordered_map(node_gmsh_to_libmesh_map,node_libmesh_to_gmsh_map);
//		carl::invert_index_unordered_map(element_gmsh_to_libmesh_map,element_libmesh_to_gmsh_map);
////
////		// DEBUG Small bcast info test
////		MeshComm.barrier();
////		std::cout << " -> " << rank << " "
////							<< node_gmsh_to_libmesh_map.size() << " / "
////							<< node_libmesh_to_gmsh_map.size() << " --- "
////							<< element_gmsh_to_libmesh_map.size() << " / "
////							<< element_libmesh_to_gmsh_map.size() <<  std::endl;
////
////		MeshComm.barrier();
////		if(rank == 0)
////		{
////			std::cout << std::endl;
////		}
////
////		MeshComm.barrier();
////		std::cout << " -> " << rank << " "
////							<< "7" << " "
////							<< element_gmsh_to_libmesh_map[7] << " "
////							<< element_libmesh_to_gmsh_map[element_gmsh_to_libmesh_map[7]]
////							<< std::endl;
////
////		MeshComm.barrier();
////		if(rank == 0)
////		{
////			std::cout << std::endl;
////		}
////		MeshComm.barrier();
//	//}
//	//else
//	//{
//	//	libmesh_error_msg("Error: unknown filetype");
//	//}
//};
//
//void carl::build_mesh_map_Gmsh(const std::string &filename, std::unordered_map<int,int> &node_map, std::unordered_map<int,int> &element_map)
//{
//	//if(!(boost::filesystem::path(filename).extension().string().compare(".msh")==0))
//	//{
//	//	libmesh_error_msg("Error: expected Gmsh filetype");
//	//}
//
//	// Open file
//	std::ifstream dataF(filename);
//	if(!dataF.good())
//	{
//		libmesh_error_msg("Error: bad filestream");
//	}
//
//	// Buffer string
//	std::string bufferLine;
//	std::stringstream	dataBuffer;
//
//	int nbOfNodes = -1;
//	int nbOfElements = -1;
//
//	int gmshNodeIndex = -1;
//	int gmshElemIndex = -1;
//
//	bool hasNodes = false;
//	bool hasElements = false;
//
//	// Read info until the file ends
//	while(std::getline(dataF,bufferLine))
//	{
//		/*
//		 * 		As of the current Gmsh version (2.11, filetype v. 2.2), each
//		 * 	file has only one $Nodes and one $Elements sections
//		 */
//		if(bufferLine.compare("$Nodes")==0)
//		{
//			// Read node indexes
//			hasNodes = true;
//			dataF >> nbOfNodes;
//
//			node_map.reserve(nbOfNodes);
//
//			// Line structure
//			// [index] [X] [Y] [Z]
//			std::getline(dataF,bufferLine);
//			for(unsigned int iii = 0; iii < nbOfNodes; ++iii)
//			{
//				std::getline(dataF,bufferLine);
//				dataBuffer.str("");
//				dataBuffer.clear();
//				dataBuffer << bufferLine;
//
//				dataBuffer >> gmshNodeIndex;
//
//				// Add point to map
//				node_map[gmshNodeIndex] = iii;
//			}
//		}
//
//		if(bufferLine.compare("$Elements")==0)
//		{
//			// Read element indexes
//			hasElements = true;
//			dataF >> nbOfElements;
//
//			element_map.reserve(nbOfElements);
//
//			// Line structure
//			// [index] [X] [Y] [Z]
//			std::getline(dataF,bufferLine);
//			for(unsigned int iii = 0; iii < nbOfElements; ++iii)
//			{
//				std::getline(dataF,bufferLine);
//				dataBuffer.str("");
//				dataBuffer.clear();
//				dataBuffer << bufferLine;
//
//				dataBuffer 	>> gmshElemIndex;
//
//				element_map[gmshElemIndex] = iii;
//			}
//		}
//
//		if(hasNodes && hasElements)
//		{
//			break;
//		}
//	}
//}
//
//void carl::build_mesh_map_Gmsh(const std::string &filename,
//		std::unordered_map<int,int> &node_gmsh_to_libmesh_map,
//		std::unordered_map<int,int> &node_libmesh_to_gmsh_map,
//		std::unordered_map<int,int> &element_gmsh_to_libmesh_map,
//		std::unordered_map<int,int> &element_libmesh_to_gmsh_map)
//{
//	//if(!(boost::filesystem::path(filename).extension().string().compare(".msh")==0))
//	//{
//	//	libmesh_error_msg("Error: expected Gmsh filetype");
//	//}
//
//	// Open file
//	std::ifstream dataF(filename);
//	if(!dataF.good())
//	{
//		libmesh_error_msg("Error: bad filestream");
//	}
//
//	// Buffer string
//	std::string bufferLine;
//	std::stringstream	dataBuffer;
//
//	int nbOfNodes = -1;
//	int nbOfElements = -1;
//
//	int gmshNodeIndex = -1;
//	int gmshElemIndex = -1;
//
//	bool hasNodes = false;
//	bool hasElements = false;
//
//	// Read info until the file ends
//	while(std::getline(dataF,bufferLine))
//	{
//		/*
//		 * 		As of the current Gmsh version (2.11, filetype v. 2.2), each
//		 * 	file has only one $Nodes and one $Elements sections
//		 */
//		if(bufferLine.compare("$Nodes")==0)
//		{
//			// Read node indexes
//			hasNodes = true;
//			dataF >> nbOfNodes;
//
//			node_gmsh_to_libmesh_map.reserve(nbOfNodes);
//			node_libmesh_to_gmsh_map.reserve(nbOfNodes);
//
//			// Line structure
//			// [index] [X] [Y] [Z]
//			std::getline(dataF,bufferLine);
//			for(unsigned int iii = 0; iii < nbOfNodes; ++iii)
//			{
//				std::getline(dataF,bufferLine);
//				dataBuffer.str("");
//				dataBuffer.clear();
//				dataBuffer << bufferLine;
//
//				dataBuffer >> gmshNodeIndex;
//
//				// Add point to map
//				node_gmsh_to_libmesh_map[gmshNodeIndex] = iii;
//				node_libmesh_to_gmsh_map[iii] = gmshNodeIndex;
//			}
//		}
//
//		if(bufferLine.compare("$Elements")==0)
//		{
//			// Read element indexes
//			hasElements = true;
//			dataF >> nbOfElements;
//
//			element_gmsh_to_libmesh_map.reserve(nbOfElements);
//			element_libmesh_to_gmsh_map.reserve(nbOfElements);
//
//			// Line structure
//			// [index] [X] [Y] [Z]
//			std::getline(dataF,bufferLine);
//			for(unsigned int iii = 0; iii < nbOfElements; ++iii)
//			{
//				std::getline(dataF,bufferLine);
//				dataBuffer.str("");
//				dataBuffer.clear();
//				dataBuffer << bufferLine;
//
//				dataBuffer 	>> gmshElemIndex;
//
//				element_gmsh_to_libmesh_map[gmshElemIndex] = iii;
//				element_libmesh_to_gmsh_map[iii] = gmshElemIndex;
//			}
//		}
//
//		if(hasNodes && hasElements)
//		{
//			break;
//		}
//	}
//}

//void carl::set_mesh_Gmsh(	libMesh::Mesh& mesh, const std::string& mesh_file,
//					std::unordered_map<int,int>& mesh_NodeMap, std::unordered_map<int,int>& mesh_ElemMap)
//{
//	libMesh::GmshIO meshBuffer(mesh);
//	meshBuffer.read(mesh_file);
//	mesh.prepare_for_use();
//	create_mesh_map(mesh_file,mesh_NodeMap,mesh_ElemMap);
//};
//
//template<typename libMeshMesh>
//void carl::set_mesh_Gmsh(	libMeshMesh& mesh, const std::string& mesh_file)
//{
//	libMesh::GmshIO meshBuffer(mesh);
//	meshBuffer.read(mesh_file);
//	mesh.prepare_for_use();
//};
//
//template void carl::set_mesh_Gmsh(	libMesh::Mesh& mesh, const std::string& mesh_file);
//template void carl::set_mesh_Gmsh(	libMesh::ParallelMesh& mesh, const std::string& mesh_file);
//template void carl::set_mesh_Gmsh(	libMesh::SerialMesh& mesh, const std::string& mesh_file);
//
//void carl::set_mesh_Gmsh(
//		libMesh::Mesh& mesh,
//		const std::string& mesh_file,
//		std::unordered_map<int,int> &node_gmsh_to_libmesh_map,
//		std::unordered_map<int,int> &node_libmesh_to_gmsh_map,
//		std::unordered_map<int,int> &element_gmsh_to_libmesh_map,
//		std::unordered_map<int,int> &element_libmesh_to_gmsh_map
//		)
//{
//	libMesh::GmshIO meshBuffer(mesh);
//	meshBuffer.read(mesh_file);
//	mesh.prepare_for_use();
//	create_mesh_map(mesh_file,node_gmsh_to_libmesh_map,node_libmesh_to_gmsh_map,
//			element_gmsh_to_libmesh_map,element_libmesh_to_gmsh_map);
//};

// SERIAL intersection tables generation
void carl::generate_intersection_tables_partial(	std::string& intersection_table_restrict_B_Filename,
											std::string& intersection_table_I_Filename,
											std::unordered_map<int,int>& mesh_restrict_ElemMap,
											std::unordered_map<int,int>& mesh_micro_ElemMap,
											std::unordered_map<int,int>& mesh_inter_ElemMap,
											std::vector<std::pair<int,int> >& intersection_table_restrict_B,
											std::unordered_multimap<int,int>& intersection_table_I
											)
{
	std::ifstream table_restrict_B_file(intersection_table_restrict_B_Filename);
	std::ifstream table_I_file(intersection_table_I_Filename);

	int nbOfIntersections_restrict_B = -1;
	int nbOfIntersectionsI  = -1;
	int nbOfTotalTetrasI  = -1;
	int dummyInt = -1;
	int nbOfTetras = -1;
	int tetraIdx = -1;

	int extIdxA = -1;
	int extIdxB = -1;
	int extIdxI = -1;

	table_restrict_B_file >> nbOfIntersections_restrict_B;
	table_I_file >> nbOfIntersectionsI >> nbOfTotalTetrasI;

	homemade_assert_msg(nbOfIntersections_restrict_B == nbOfIntersectionsI, "Incompatible intersection table files!");

	intersection_table_restrict_B.resize(nbOfIntersections_restrict_B);
	intersection_table_I.reserve(2*nbOfTotalTetrasI);

	for(int iii = 0; iii < nbOfIntersections_restrict_B; ++iii)
	{
		table_restrict_B_file >> dummyInt >> extIdxA >> extIdxB;
		intersection_table_restrict_B[iii].first = mesh_restrict_ElemMap[extIdxA];
		intersection_table_restrict_B[iii].second = mesh_micro_ElemMap[extIdxB];

		table_I_file >> dummyInt >> nbOfTetras;
		for(int jjj = 0; jjj < nbOfTetras; ++jjj)
		{
			table_I_file >> extIdxI;
			tetraIdx = mesh_inter_ElemMap[extIdxI];
			intersection_table_I.emplace(dummyInt,tetraIdx);
		}
	}

	table_restrict_B_file.close();
	table_I_file.close();
};

void carl::generate_intersection_tables_full(		std::string& equivalence_table_restrict_A_Filename,
											std::string& intersection_table_restrict_B_Filename,
											std::string& intersection_table_I_Filename,
											std::unordered_map<int,int>& mesh_restrict_ElemMap,
											std::unordered_map<int,int>& mesh_micro_ElemMap,
											std::unordered_map<int,int>& mesh_BIG_ElemMap,
											std::unordered_map<int,int>& mesh_inter_ElemMap,
											std::unordered_map<int,int>& equivalence_table_restrict_A,
											std::vector<std::pair<int,int> >& intersection_table_restrict_B,
											std::unordered_multimap<int,int>& intersection_table_I
											)
{
	std::ifstream table_restrict_A_file(equivalence_table_restrict_A_Filename);
	std::ifstream table_restrict_B_file(intersection_table_restrict_B_Filename);
	std::ifstream table_I_file(intersection_table_I_Filename);

	int nbOfEquivalences_restrict_A = -1;
	int nbOfIntersections_restrict_B = -1;
	int nbOfIntersectionsI  = -1;
	int nbOfTotalTetrasI  = -1;
	int dummyInt = -1;
	int nbOfTetras = -1;
	int tetraIdx = -1;

	int extIdxA = -1;
	int extIdxB = -1;
	int extIdxI = -1;
	int extIdxR = -1;

	int idxRestrict = -1;
	int idxA = -1;

	table_restrict_A_file >> nbOfEquivalences_restrict_A;
	std::unordered_map<int,int> temp_equivalence_table_A_restrict(nbOfEquivalences_restrict_A);
	equivalence_table_restrict_A.reserve(nbOfEquivalences_restrict_A);

	for(int iii = 0; iii < nbOfEquivalences_restrict_A; ++iii)
	{
		table_restrict_A_file >> extIdxR >> extIdxA;

		idxA = mesh_BIG_ElemMap[extIdxA];
		idxRestrict = mesh_restrict_ElemMap[extIdxR];

		temp_equivalence_table_A_restrict[idxA] = idxRestrict;
		equivalence_table_restrict_A[idxRestrict] = idxA;
	}

	table_restrict_B_file >> nbOfIntersections_restrict_B;
	table_I_file >> nbOfIntersectionsI >> nbOfTotalTetrasI;

	homemade_assert_msg(nbOfIntersections_restrict_B == nbOfIntersectionsI, "Incompatible intersection table files!");

	intersection_table_restrict_B.resize(nbOfIntersections_restrict_B);
	intersection_table_I.reserve(nbOfTotalTetrasI);

	for(int iii = 0; iii < nbOfIntersections_restrict_B; ++iii)
	{
		table_restrict_B_file >> dummyInt >> extIdxA >> extIdxB;

		idxA = mesh_BIG_ElemMap[extIdxA];
		intersection_table_restrict_B[iii].first = temp_equivalence_table_A_restrict[idxA];
		intersection_table_restrict_B[iii].second = mesh_micro_ElemMap[extIdxB];

		table_I_file >> dummyInt >> nbOfTetras;
		for(int jjj = 0; jjj < nbOfTetras; ++jjj)
		{
			table_I_file >> extIdxI;
			tetraIdx = mesh_inter_ElemMap[extIdxI];
			intersection_table_I.emplace(dummyInt,tetraIdx);
		}
	}

	table_restrict_A_file.close();
	table_restrict_B_file.close();
	table_I_file.close();
};

// PARALLEL intersection tables generation
void carl::build_intersection_and_restriction_tables(
		const libMesh::Parallel::Communicator& WorldComm,
		const std::string& intersection_full_table_Filename,
		const std::string& equivalence_table_A_Filename,
		const std::string& equivalence_table_B_Filename,
		std::vector<carl::IntersectionData>& intersection_full_table,
		std::unordered_map<int,int>& equivalence_table_A_to_R_A,
		std::unordered_map<int,int>& equivalence_table_B_to_R_B,
		std::unordered_map<int,int>& equivalence_table_R_A_to_A,
		std::unordered_map<int,int>& equivalence_table_R_B_to_B
		)
{
	int rank = WorldComm.rank();

	// 	While the equivalence tables are saved as unordered maps, it's easier to
	// save them as vectors at first, broadcast them, and them reconvert to maps
	std::vector<int> dummy_equivalence_table_A;
	std::vector<int> dummy_equivalence_table_B;
	std::vector<int> dummy_intersection_full_table;

	int nbOfIntersections = -1;
	int nbOfInterElems = -1;
	int nbOfRestricted_A_Elems = -1;
	int nbOfRestricted_B_Elems = -1;

	// Do the file reading job on proc 0
	if(rank == 0)
	{
		std::ifstream intersection_full_file(intersection_full_table_Filename);

		intersection_full_file >> nbOfIntersections >> nbOfInterElems;
		dummy_intersection_full_table.resize(4*nbOfInterElems);

		for(int iii = 0; iii < nbOfInterElems; ++iii)
		{
			intersection_full_file
				>> dummy_intersection_full_table[4*iii]
				>> dummy_intersection_full_table[4*iii + 1]
				>> dummy_intersection_full_table[4*iii + 2]
				>> dummy_intersection_full_table[4*iii + 3];
		}
		intersection_full_file.close();

		std::ifstream equivalence_A_file(equivalence_table_A_Filename);

		equivalence_A_file >> nbOfRestricted_A_Elems;
		dummy_equivalence_table_A.resize(2*nbOfRestricted_A_Elems);

		for(int iii = 0; iii < nbOfRestricted_A_Elems; ++iii)
		{
			equivalence_A_file	>> dummy_equivalence_table_A[2*iii]
								>> dummy_equivalence_table_A[2*iii + 1];
		}
		equivalence_A_file.close();

		std::ifstream equivalence_B_file(equivalence_table_B_Filename);

		equivalence_B_file >> nbOfRestricted_B_Elems;
		dummy_equivalence_table_B.resize(2*nbOfRestricted_B_Elems);

		for(int iii = 0; iii < nbOfRestricted_B_Elems; ++iii)
		{
			equivalence_B_file	>> dummy_equivalence_table_B[2*iii]
								>> dummy_equivalence_table_B[2*iii + 1];
		}
		equivalence_B_file.close();
	}

	// Broadcast the sizes and resize on other procs
	WorldComm.broadcast(nbOfIntersections);
	WorldComm.broadcast(nbOfInterElems);
	WorldComm.broadcast(nbOfRestricted_A_Elems);
	WorldComm.broadcast(nbOfRestricted_B_Elems);

	if(rank != 0)
	{
		dummy_intersection_full_table.resize(4*nbOfInterElems);
		dummy_equivalence_table_A.resize(2*nbOfRestricted_A_Elems);
		dummy_equivalence_table_B.resize(2*nbOfRestricted_B_Elems);
	}

	WorldComm.barrier();

	WorldComm.broadcast(dummy_intersection_full_table);
	WorldComm.broadcast(dummy_equivalence_table_A);
	WorldComm.broadcast(dummy_equivalence_table_B);

	// Convert back to unoredered maps and pair vectors
	intersection_full_table.resize(nbOfInterElems);
	equivalence_table_A_to_R_A.reserve(nbOfRestricted_A_Elems);
	equivalence_table_B_to_R_B.reserve(nbOfRestricted_B_Elems);

	equivalence_table_R_A_to_A.reserve(nbOfRestricted_A_Elems);
	equivalence_table_R_B_to_B.reserve(nbOfRestricted_B_Elems);

	for(int iii = 0; iii < nbOfInterElems; ++iii)
	{
		intersection_full_table[iii].InterMeshIdx
									= dummy_intersection_full_table[4*iii];
		intersection_full_table[iii].AMeshIdx
									= dummy_intersection_full_table[4*iii + 1];
		intersection_full_table[iii].BMeshIdx
									= dummy_intersection_full_table[4*iii + 2];
		intersection_full_table[iii].IntersectionID
									= dummy_intersection_full_table[4*iii + 3];
	}

	for(int iii = 0; iii < nbOfRestricted_A_Elems; ++iii)
	{
		equivalence_table_R_A_to_A[dummy_equivalence_table_A[2*iii]] =
				dummy_equivalence_table_A[2*iii + 1];
		equivalence_table_A_to_R_A[dummy_equivalence_table_A[2*iii + 1]] =
				dummy_equivalence_table_A[2*iii];
	}

	for(int iii = 0; iii < nbOfRestricted_B_Elems; ++iii)
	{
		equivalence_table_R_B_to_B[dummy_equivalence_table_B[2*iii]] =
				dummy_equivalence_table_B[2*iii + 1];
		equivalence_table_B_to_R_B[dummy_equivalence_table_B[2*iii + 1]] =
				dummy_equivalence_table_B[2*iii];
	}

//	// DEBUG Small bcast info test
//	std::cout << " -> " << rank << " "
//						<< intersection_full_table.size() << " / "
//						<< nbOfInterElems << " --- "
//						<< dummy_equivalence_table_A.size()/2 << " / "
//						<< nbOfRestricted_A_Elems << " --- "
//						<< dummy_equivalence_table_B.size()/2 << " / "
//						<< nbOfRestricted_B_Elems << std::endl;
//
//	WorldComm.barrier();
//	if(rank == 0)
//	{
//		std::cout << std::endl;
//	}
//	WorldComm.barrier();
//	std::cout << " -> " << rank << " "
//						<< intersection_full_table[4].InterMeshIdx << " "
//						<< intersection_full_table[4].AMeshIdx << " "
//						<< intersection_full_table[4].BMeshIdx << " "
//						<< intersection_full_table[4].IntersectionID << " "
//						<< std::endl;
//
//	WorldComm.barrier();
//	if(rank == 0)
//	{
//		std::cout << std::endl;
//	}
//	WorldComm.barrier();
//	std::cout << " -> " << rank << " "
//						<< "4" << " "
//						<< equivalence_table_A_to_R_A[4] << " "
//						<< equivalence_table_R_A_to_A[equivalence_table_A_to_R_A[4]]
//						<< std::endl;
//
//	WorldComm.barrier();
//	if(rank == 0)
//	{
//		std::cout << std::endl;
//	}
//	WorldComm.barrier();
//	std::cout << " -> " << rank << " "
//						<< "4" << " "
//						<< equivalence_table_B_to_R_B[4] << " "
//						<< equivalence_table_R_B_to_B[equivalence_table_B_to_R_B[4]]
//						<< std::endl;
};

void carl::set_equivalence_tables(
		const libMesh::Parallel::Communicator& WorldComm,
		const std::string& equivalence_table_A_Filename,
		const std::string& equivalence_table_B_Filename,

		std::unordered_map<int,int>& equivalence_table_A_to_R_A,
		std::unordered_map<int,int>& equivalence_table_B_to_R_B,
		std::unordered_map<int,int>& equivalence_table_R_A_to_A,
		std::unordered_map<int,int>& equivalence_table_R_B_to_B )
{
	int rank = WorldComm.rank();

	// 	While the equivalence tables are saved as unordered maps, it's easier to
	// save them as vectors at first, broadcast them, and them reconvert to maps
	std::vector<int> dummy_equivalence_table_A;
	std::vector<int> dummy_equivalence_table_B;

	int nbOfRestricted_A_Elems = -1;
	int nbOfRestricted_B_Elems = -1;

	// Read files with proc 0
	if(rank == 0)
	{
		int temp_RX = -1;
		int temp_X = -1;

		std::ifstream equivalence_A_file(equivalence_table_A_Filename);

		equivalence_A_file >> nbOfRestricted_A_Elems;
		dummy_equivalence_table_A.resize(2*nbOfRestricted_A_Elems);

		for(int iii = 0; iii < nbOfRestricted_A_Elems; ++iii)
		{
			equivalence_A_file	>> temp_RX
								>> temp_X;

			dummy_equivalence_table_A[2*iii] = temp_RX;
			dummy_equivalence_table_A[2*iii + 1] = temp_X;
		}
		equivalence_A_file.close();

		std::ifstream equivalence_B_file(equivalence_table_B_Filename);

		equivalence_B_file >> nbOfRestricted_B_Elems;
		dummy_equivalence_table_B.resize(2*nbOfRestricted_B_Elems);

		for(int iii = 0; iii < nbOfRestricted_B_Elems; ++iii)
		{
			equivalence_B_file	>> temp_RX
								>> temp_X;

			dummy_equivalence_table_B[2*iii] = temp_RX;
			dummy_equivalence_table_B[2*iii + 1] = temp_X;
		}
		equivalence_B_file.close();
	}

	// Broadcast the sizes and resize on other procs
	WorldComm.broadcast(nbOfRestricted_A_Elems);
	WorldComm.broadcast(nbOfRestricted_B_Elems);

	if(rank != 0)
	{
		dummy_equivalence_table_A.resize(2*nbOfRestricted_A_Elems);
		dummy_equivalence_table_B.resize(2*nbOfRestricted_B_Elems);
	}

	WorldComm.barrier();

	WorldComm.broadcast(dummy_equivalence_table_A);
	WorldComm.broadcast(dummy_equivalence_table_B);

	// Convert back to unoredered maps
	equivalence_table_A_to_R_A.reserve(nbOfRestricted_A_Elems);
	equivalence_table_B_to_R_B.reserve(nbOfRestricted_B_Elems);

	equivalence_table_R_A_to_A.reserve(nbOfRestricted_A_Elems);
	equivalence_table_R_B_to_B.reserve(nbOfRestricted_B_Elems);

	for(int iii = 0; iii < nbOfRestricted_A_Elems; ++iii)
	{
		equivalence_table_R_A_to_A[dummy_equivalence_table_A[2*iii]] =
				dummy_equivalence_table_A[2*iii + 1];
		equivalence_table_A_to_R_A[dummy_equivalence_table_A[2*iii + 1]] =
				dummy_equivalence_table_A[2*iii];
	}

	for(int iii = 0; iii < nbOfRestricted_B_Elems; ++iii)
	{
		equivalence_table_R_B_to_B[dummy_equivalence_table_B[2*iii]] =
				dummy_equivalence_table_B[2*iii + 1];
		equivalence_table_B_to_R_B[dummy_equivalence_table_B[2*iii + 1]] =
				dummy_equivalence_table_B[2*iii];
	}

//	// DEBUG Small bcast info test
//	std::cout << " -> " << rank << " "
//						<< dummy_equivalence_table_A.size()/2 << " / "
//						<< nbOfRestricted_A_Elems << " --- "
//						<< dummy_equivalence_table_B.size()/2 << " / "
//						<< nbOfRestricted_B_Elems << std::endl;
//
//	WorldComm.barrier();
//	if(rank == 0)
//	{
//		std::cout << std::endl;
//	}
//
//	WorldComm.barrier();
//	std::cout << " -> " << rank << " "
//						<< "7" << " "
//						<< equivalence_table_R_A_to_A[7] << " "
//						<< equivalence_table_A_to_R_A[equivalence_table_R_A_to_A[7]]
//						<< std::endl;
//
//	WorldComm.barrier();
//	if(rank == 0)
//	{
//		std::cout << std::endl;
//	}
//	WorldComm.barrier();
//	std::cout << " -> " << rank << " "
//						<< "4" << " "
//						<< equivalence_table_R_B_to_B[4] << " "
//						<< equivalence_table_B_to_R_B[equivalence_table_R_B_to_B[4]]
//						<< std::endl;
//	WorldComm.barrier();
//	if(rank == 0)
//	{
//		std::cout << std::endl;
//	}
//	WorldComm.barrier();
}

void carl::set_restricted_intersection_pairs_table(
		const std::unordered_map<int,std::pair<int,int> >&  full_intersection_pairs_map,
		const std::unordered_map<int,int>& equivalence_table_A_to_R_A,
		const std::unordered_map<int,int>& equivalence_table_B_to_R_B,
		std::unordered_map<int,std::pair<int,int> >& full_intersection_restricted_pairs_map)
{
	// "Resize" the final output
	full_intersection_restricted_pairs_map.reserve(equivalence_table_A_to_R_A.size());

	int interID = -1;
	int idxA = -1;
	int idxB = -1;

	int idxRA = -1;
	int idxRB = -1;

	std::unordered_map<int,std::pair<int,int> >::const_iterator mapIt =
			full_intersection_pairs_map.begin();
	std::unordered_map<int,std::pair<int,int> >::const_iterator end_mapIt =
			full_intersection_pairs_map.end();

	for( ; mapIt != end_mapIt; ++mapIt)
	{
		interID	= mapIt->first;
		idxA	= mapIt->second.first;
		idxB	= mapIt->second.second;

		idxRA	= equivalence_table_A_to_R_A.at(idxA);
		idxRB	= equivalence_table_B_to_R_B.at(idxB);

		full_intersection_restricted_pairs_map[interID] =
				std::pair<int,int>(idxRA,idxRB);
	}
}

void carl::set_full_intersection_tables(
		const libMesh::Parallel::Communicator& WorldComm,
		const std::string& intersection_full_table_Filename,

		std::unordered_map<int,std::pair<int,int> >& full_intersection_pairs_map,
		std::unordered_map<int,int>& full_intersection_meshI_to_inter_map)
{
	int rank = WorldComm.rank();

	// 	While the equivalence tables are saved as unordered maps, it's easier to
	// save them as vectors at first, broadcast them, and them reconvert to maps
	std::vector<int> dummy_intersections_IDs;
	std::vector<int> dummy_intersection_pairs_table;
	std::vector<int> dummy_all_intersections_table;
	std::vector<int> dummy_intersections_sizes_table;

	int nbOfIntersections = -1;
	int nbOfInterElems = -1;

	// Declare a few auxiliary variables
	int temp_interID = -1;
	int temp_idxA = -1;
	int temp_idxB = -1;
	int temp_idxI = -1;
	int temp_nbOfInter = -1;

	int interIdx = 0;

	// Do the file reading job on proc 0
	if(rank == 0)
	{
		std::ifstream intersection_full_file(intersection_full_table_Filename);

		intersection_full_file >> nbOfIntersections >> nbOfInterElems;

		dummy_intersections_IDs.resize(nbOfIntersections);
		dummy_intersection_pairs_table.resize(2*nbOfIntersections);
		dummy_all_intersections_table.resize(nbOfInterElems);
		dummy_intersections_sizes_table.resize(nbOfIntersections);

		for(int iii = 0; iii < nbOfIntersections; ++iii)
		{
			// For each line, read ...

			// ... the inter ID, A and B's elements, and the number of I's elements ...
			intersection_full_file 	>> temp_interID
									>> temp_idxA >> temp_idxB
									>> temp_nbOfInter;

			dummy_intersections_IDs[iii]			  = temp_interID;
			dummy_intersection_pairs_table[2*iii]     = temp_idxA;
			dummy_intersection_pairs_table[2*iii + 1] = temp_idxB;
			dummy_intersections_sizes_table[iii]      = temp_nbOfInter;

			// ... and all of I's elements
			for(int jjj = 0; jjj < temp_nbOfInter; ++jjj)
			{
				intersection_full_file >> temp_idxI;
				dummy_all_intersections_table[interIdx] = temp_idxI;

				++interIdx;
			}
		}
		intersection_full_file.close();
	}

	// Broadcast the sizes and resize on other procs
	WorldComm.broadcast(nbOfIntersections);
	WorldComm.broadcast(nbOfInterElems);

	if(rank != 0)
	{
		dummy_intersections_IDs.resize(nbOfIntersections);
		dummy_intersection_pairs_table.resize(2*nbOfIntersections);
		dummy_all_intersections_table.resize(nbOfInterElems);
		dummy_intersections_sizes_table.resize(nbOfIntersections);
	}

	WorldComm.barrier();

	WorldComm.broadcast(dummy_intersections_IDs);
	WorldComm.broadcast(dummy_intersection_pairs_table);
	WorldComm.broadcast(dummy_all_intersections_table);
	WorldComm.broadcast(dummy_intersections_sizes_table);

	full_intersection_pairs_map.reserve(nbOfIntersections);
	full_intersection_meshI_to_inter_map.reserve(nbOfInterElems);

	interIdx = 0;

	// Convert the vectors to maps
	for(int iii = 0; iii < nbOfIntersections; ++iii)
	{
		temp_interID 	= dummy_intersections_IDs[iii];
		temp_idxA		= dummy_intersection_pairs_table[2*iii];
		temp_idxB		= dummy_intersection_pairs_table[2*iii + 1];

		full_intersection_pairs_map[temp_interID] = std::pair<int,int>(temp_idxA,temp_idxB);

		temp_nbOfInter  = dummy_intersections_sizes_table[iii];
		for(int jjj = 0; jjj < temp_nbOfInter; ++jjj)
		{
			temp_idxI	= dummy_all_intersections_table[interIdx];
			full_intersection_meshI_to_inter_map[temp_idxI] = temp_interID;
			++interIdx;
		}
	}

//	// DEBUG Small bcast info test
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD," -> %d %d / %d\n",rank,full_intersection_meshI_to_inter_map.size(),nbOfInterElems);
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//	WorldComm.barrier();
//	if(rank == 0)
//	{
//		std::cout << std::endl;
//	}
//	WorldComm.barrier();
//
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD," -> %d %d / %d\n",rank,full_intersection_pairs_map.size(),nbOfIntersections);
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//	WorldComm.barrier();
//	if(rank == 0)
//	{
//		std::cout << std::endl;
//	}
//	WorldComm.barrier();
//
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD," -> %d %d %d %d \n",
//			rank,full_intersection_meshI_to_inter_map[7],
//			full_intersection_pairs_map[full_intersection_meshI_to_inter_map[7]].first,
//			full_intersection_pairs_map[full_intersection_meshI_to_inter_map[7]].second);
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//	WorldComm.barrier();
//	if(rank == 0)
//	{
//		std::cout << std::endl;
//	}
//	WorldComm.barrier();
};

void carl::set_intersection_tables(
		const libMesh::Parallel::Communicator& WorldComm,
		const libMesh::Mesh& mesh_intersection,
		const std::string& intersection_full_table_Filename,
		const std::string& equivalence_table_A_Filename,
		const std::string& equivalence_table_B_Filename,

		const std::unordered_map<int,int>& equivalence_table_A_to_R_A,
		const std::unordered_map<int,int>& equivalence_table_B_to_R_B,

		std::unordered_map<int,std::pair<int,int> >& full_intersection_pairs_map,
		std::unordered_map<int,std::pair<int,int> >& full_intersection_restricted_pairs_map,
		std::unordered_map<int,int>& local_intersection_meshI_to_inter_map
		)
{
	int rank = WorldComm.rank();

	// 	Start by reading and broadcasting the global intersection table
	//  TODO : 	change this step after parallelizing and integrating the
	//			intersection algorithm
	std::unordered_map<int,int> full_intersection_meshI_to_inter_map;
	set_full_intersection_tables(WorldComm,intersection_full_table_Filename,
			full_intersection_pairs_map,full_intersection_meshI_to_inter_map);

	// Set up the restricted intersection pairs table
	set_restricted_intersection_pairs_table(full_intersection_pairs_map,
		equivalence_table_A_to_R_A,equivalence_table_B_to_R_B,
		full_intersection_restricted_pairs_map);

//	// DEBUG Print all this into a file
//	{
//		std::string debug_inter_table = "debug_restriction_libmesh_" +
//				std::to_string(rank) + ".dat";
//		std::ofstream debug_inter_file(debug_inter_table);
//
//		std::unordered_map<int,std::pair<int,int> >::iterator mapIt = full_intersection_restricted_pairs_map.begin();
//		std::unordered_map<int,std::pair<int,int> >::iterator end_mapIt = full_intersection_restricted_pairs_map.end();
//
//		for( ; mapIt != end_mapIt; ++mapIt)
//		{
//			debug_inter_file << mapIt->first << " " << mapIt->second.first << " " << mapIt->second.second << std::endl;
//		}
//
//		debug_inter_file.close();
//	}

	// Build in each processor its own intersection table
	libMesh::MeshBase::const_element_iterator       elemIt  = mesh_intersection.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_elemIt = mesh_intersection.active_local_elements_end();
	int local_nbOfInters = mesh_intersection.n_active_local_elem();
	local_intersection_meshI_to_inter_map.reserve(local_nbOfInters);

	// Some dummy indexes used
	int idxI_table = -1;
	int interID = -1;

	for( ; elemIt != end_elemIt; ++ elemIt)
	{
		const libMesh::Elem* elem = *elemIt;
		idxI_table 	= elem->id();
		interID 	= full_intersection_meshI_to_inter_map[idxI_table];
		local_intersection_meshI_to_inter_map[idxI_table] = interID;
	}

//	// DEBUG Print all this into a file
//	{
//		std::string debug_inter_table = "debug_local_inter_libmesh_" +
//				std::to_string(rank) + ".dat";
//		std::ofstream debug_inter_file(debug_inter_table);
//
//		std::unordered_map<int,int>::iterator mapIt = local_intersection_meshI_to_inter_map.begin();
//		std::unordered_map<int,int>::iterator end_mapIt = local_intersection_meshI_to_inter_map.end();
//
//		for( ; mapIt != end_mapIt; ++mapIt)
//		{
//			debug_inter_file << mapIt->first << " " << mapIt->second << std::endl;
//		}
//
//		debug_inter_file.close();
//	}
};
