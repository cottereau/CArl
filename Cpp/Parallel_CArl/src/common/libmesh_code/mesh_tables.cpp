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
}

void carl::set_restricted_intersection_pairs_table(
		const std::unordered_map<int,std::pair<int,int> >&  intersection_pairs_map,
		const std::unordered_map<int,int>& equivalence_table_A_to_R_A,
		const std::unordered_map<int,int>& equivalence_table_B_to_R_B,
		std::unordered_map<int,std::pair<int,int> >& intersection_restricted_pairs_map)
{
	// "Resize" the final output
	intersection_restricted_pairs_map.reserve(equivalence_table_A_to_R_A.size());

	int interID = -1;
	int idxA = -1;
	int idxB = -1;

	int idxRA = -1;
	int idxRB = -1;

	std::unordered_map<int,std::pair<int,int> >::const_iterator mapIt =
			intersection_pairs_map.begin();
	std::unordered_map<int,std::pair<int,int> >::const_iterator end_mapIt =
			intersection_pairs_map.end();

	for( ; mapIt != end_mapIt; ++mapIt)
	{
		interID	= mapIt->first;
		idxA	= mapIt->second.first;
		idxB	= mapIt->second.second;

		idxRA	= equivalence_table_A_to_R_A.at(idxA);
		idxRB	= equivalence_table_B_to_R_B.at(idxB);

		intersection_restricted_pairs_map[interID] =
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
};

void carl::read_local_intersection_tables(
		const libMesh::Parallel::Communicator& WorldComm,
		const std::string& intersection_local_table_Filename,

		std::unordered_map<int,std::pair<int,int> >& local_intersection_pairs_map,
		std::unordered_map<int,int>& local_intersection_meshI_to_inter_map)
{
	int nbOfIntersections = -1;
	int nbOfInterElems = -1;

	// Declare a few auxiliary variables
	int temp_interID = -1;
	int temp_idxA = -1;
	int temp_idxB = -1;
	int temp_idxI = -1;
	int temp_nbOfInter = -1;
	int temp_dummy = -1;

	int interIdx = 0;

	// Do the file reading
	std::ifstream intersection_full_file(intersection_local_table_Filename);

	intersection_full_file >> nbOfIntersections >> nbOfInterElems >> temp_dummy;

	for(int iii = 0; iii < nbOfIntersections; ++iii)
	{
		// For each line, read ...

		// ... the inter ID, A and B's elements, and the number of I's elements ...
		intersection_full_file 	>> temp_interID
								>> temp_idxA >> temp_idxB
								>> temp_nbOfInter;

		local_intersection_pairs_map[temp_interID] = std::pair<int,int>(temp_idxA,temp_idxB);

		// ... and all of I's elements
		for(int jjj = 0; jjj < temp_nbOfInter; ++jjj)
		{
			intersection_full_file >> temp_idxI;
			local_intersection_meshI_to_inter_map[temp_idxI] = temp_interID;
			++interIdx;
		}
	}
	intersection_full_file.close();

	WorldComm.barrier();
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
	// 	Start by reading and broadcasting the global intersection table
	std::unordered_map<int,int> full_intersection_meshI_to_inter_map;
	set_full_intersection_tables(WorldComm,intersection_full_table_Filename,
			full_intersection_pairs_map,full_intersection_meshI_to_inter_map);

	// Set up the restricted intersection pairs table
	set_restricted_intersection_pairs_table(full_intersection_pairs_map,
		equivalence_table_A_to_R_A,equivalence_table_B_to_R_B,
		full_intersection_restricted_pairs_map);

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
};

void carl::set_local_intersection_tables(
		const libMesh::Parallel::Communicator& WorldComm,
		const libMesh::Mesh& mesh_intersection,
		const std::string& intersection_local_table_Filename,
		const std::string& equivalence_table_A_Filename,
		const std::string& equivalence_table_B_Filename,

		const std::unordered_map<int,int>& equivalence_table_A_to_R_A,
		const std::unordered_map<int,int>& equivalence_table_B_to_R_B,

		std::unordered_map<int,std::pair<int,int> >& local_intersection_pairs_map,
		std::unordered_map<int,std::pair<int,int> >& local_intersection_restricted_pairs_map,
		std::unordered_map<int,int>& local_intersection_meshI_to_inter_map
		)
{
	// 	Start by reading and broadcasting the global intersection table
	std::unordered_map<int,int> full_intersection_meshI_to_inter_map;
	read_local_intersection_tables(WorldComm,intersection_local_table_Filename,
			local_intersection_pairs_map,full_intersection_meshI_to_inter_map);

	// Set up the restricted intersection pairs table
	set_restricted_intersection_pairs_table(local_intersection_pairs_map,
		equivalence_table_A_to_R_A,equivalence_table_B_to_R_B,
		local_intersection_restricted_pairs_map);

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
};

void carl::set_global_mediator_system_intersection_lists(
		const libMesh::Parallel::Communicator& WorldComm,
		const std::string& intersection_global_table_Filename,
		const std::unordered_map<int,int>& equivalence_table_system_to_mediator,
		const std::unordered_map<int,int>& equivalence_table_mediator_to_system,

		std::unordered_multimap<int,int>& inter_mediator_A,
		std::unordered_multimap<int,int>& inter_mediator_B)
{
	// First, do the reading work on proc. 0
	int rank = WorldComm.rank();
	int nodes = WorldComm.size();

	int 				temp_interID;
	std::vector<int>	temp_idxA;
	std::vector<int>	temp_idxB;

	int nbOfIntersections = -1;
	int nbOfInterElems = -1;

	std::string dummy_data;

	if(rank == 0)
	{
		std::ifstream intersection_full_file(intersection_global_table_Filename);
		intersection_full_file >> nbOfIntersections >> nbOfInterElems;

		temp_idxA.resize(nbOfIntersections);
		temp_idxB.resize(nbOfIntersections);

		for(int iii = 0; iii < nbOfIntersections; ++iii)
		{
			intersection_full_file 	>> temp_interID
									>> temp_idxA[iii] >> temp_idxB[iii];
			intersection_full_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}
		intersection_full_file.close();
	}

	if(nodes > 1)
	{
		WorldComm.broadcast(nbOfIntersections);

		if(rank != 0)
		{
			temp_idxA.resize(nbOfIntersections);
			temp_idxB.resize(nbOfIntersections);
		}
		WorldComm.broadcast(temp_idxA);
		WorldComm.broadcast(temp_idxB);
	}

	inter_mediator_A.reserve(nbOfIntersections);
	inter_mediator_B.reserve(nbOfIntersections);

	int mediator_idx;
	for(int iii = 0; iii < nbOfIntersections; ++iii)
	{
		mediator_idx = equivalence_table_system_to_mediator.at(temp_idxA[iii]);
		inter_mediator_B.emplace(std::make_pair(mediator_idx, temp_idxB[iii]));
	}

	int system_idx;
	for(unsigned int iii = 0; iii < equivalence_table_system_to_mediator.size(); ++iii)
	{
		system_idx = equivalence_table_mediator_to_system.at(iii);
		inter_mediator_A.emplace(std::make_pair(iii, system_idx));
	}
};
