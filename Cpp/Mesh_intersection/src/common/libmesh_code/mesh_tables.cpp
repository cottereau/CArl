#include "mesh_tables.h"

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

void carl::create_mesh_map(const std::string &filename, std::unordered_map<int,int> &node_map, std::unordered_map<int,int> &element_map)
{
	/*
	 *  libmesh : nodes and elements start at 0, and are continuous
	 *  GMSH    : nodes and elements can be discontinuous, libmesh builds a map
	 *  Abaqus  : not sure, but seems similar to GMSH. libmesh builds a map, but
	 *            differently from the method used for GMSH
	 */

	if(boost::filesystem::path(filename).extension().string().compare(".msh")==0)
	{
		build_mesh_map_Gmsh(filename,node_map,element_map);
	}
	else
	{
		libmesh_error_msg("Error: unknown filetype");
	}
}

void carl::build_mesh_map_Gmsh(const std::string &filename, std::unordered_map<int,int> &node_map, std::unordered_map<int,int> &element_map)
{
	if(!(boost::filesystem::path(filename).extension().string().compare(".msh")==0))
	{
		libmesh_error_msg("Error: expected Gmsh filetype");
	}

	// Open file
	std::ifstream dataF(filename);
	if(!dataF.good())
	{
		libmesh_error_msg("Error: bad filestream");
	}

	// Buffer string
	std::string bufferLine;
	std::stringstream	dataBuffer;

	int nbOfNodes = -1;
	int nbOfElements = -1;

	int gmshNodeIndex = -1;
	int gmshElemIndex = -1;

	bool hasNodes = false;
	bool hasElements = false;

	// Read info until the file ends
	while(std::getline(dataF,bufferLine))
	{
		/*
		 * 		As of the current Gmsh version (2.11, filetype v. 2.2), each
		 * 	file has only one $Nodes and one $Elements sections
		 */
		if(bufferLine.compare("$Nodes")==0)
		{
			// Read node indexes
			hasNodes = true;
			dataF >> nbOfNodes;

			node_map.reserve(nbOfNodes);

			// Line structure
			// [index] [X] [Y] [Z]
			std::getline(dataF,bufferLine);
			for(unsigned int iii = 0; iii < nbOfNodes; ++iii)
			{
				std::getline(dataF,bufferLine);
				dataBuffer.str("");
				dataBuffer.clear();
				dataBuffer << bufferLine;

				dataBuffer >> gmshNodeIndex;

				// Add point to map
				node_map[gmshNodeIndex] = iii;
			}
		}

		if(bufferLine.compare("$Elements")==0)
		{
			// Read element indexes
			hasElements = true;
			dataF >> nbOfElements;

			element_map.reserve(nbOfElements);

			// Line structure
			// [index] [X] [Y] [Z]
			std::getline(dataF,bufferLine);
			for(unsigned int iii = 0; iii < nbOfElements; ++iii)
			{
				std::getline(dataF,bufferLine);
				dataBuffer.str("");
				dataBuffer.clear();
				dataBuffer << bufferLine;

				dataBuffer 	>> gmshElemIndex;

				element_map[gmshElemIndex] = iii;
			}
		}

		if(hasNodes && hasElements)
		{
			break;
		}
	}
}

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

void carl::set_mesh_Gmsh(	libMesh::Mesh& mesh, const std::string& mesh_file,
					std::unordered_map<int,int>& mesh_NodeMap, std::unordered_map<int,int>& mesh_ElemMap)
{
	libMesh::GmshIO meshBuffer(mesh);
	meshBuffer.read(mesh_file);
	mesh.prepare_for_use();
	create_mesh_map(mesh_file,mesh_NodeMap,mesh_ElemMap);
};

void carl::set_mesh_Gmsh(	libMesh::Mesh& mesh, const std::string& mesh_file)
{
	libMesh::GmshIO meshBuffer(mesh);
	meshBuffer.read(mesh_file);
	mesh.prepare_for_use();
};
