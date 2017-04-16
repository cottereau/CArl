/*
 * common_functions.h
 *
 *  Created on: Jan 26, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_COMMON_FUNCTIONS_H_
#define COMMON_COMMON_FUNCTIONS_H_

#include "common_header_libmesh.h"
#include "common_header.h"

int kronecker_delta(unsigned int i,
				   unsigned int j);

void clear_line();

/*
 *  A little structure used to save the intersection data
 */

namespace carl
{

void invert_index_unordered_map(
		const std::unordered_map<int,int>& input_map,
		std::unordered_map<int,int>& output_map);

template<typename T>
void jump_lines(T& filestream, unsigned int numberOfLines = 1)
{
	std::string dummy;
	for(unsigned int iii = 0; iii < numberOfLines; ++iii)
		std::getline(filestream,dummy);
};

int voigt_index_converter(int aaa, int bbb);

struct IntersectionData
{
	int InterMeshIdx;
	int AMeshIdx;
	int BMeshIdx;
	int IntersectionID;
};

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

void print_stats_to_file(std::vector<double>& vec_data, const std::string filename);

// Do libMesh's system initialization without the matrix and vector memory allocations
template<typename Sys>
void reduced_system_init(Sys& system_input)
{
	libMesh::DofMap& system_dof_map = system_input.get_dof_map();
	libMesh::MeshBase& system_mesh = system_input.get_mesh();
	
	unsigned int nb_of_variable_groups = system_input.n_variable_groups();
	for (unsigned int vg=0; vg<nb_of_variable_groups; vg++)
	{
		system_dof_map.add_variable_group(system_input.variable_group(vg));
	}

	system_dof_map.distribute_dofs(system_mesh);
	system_input.reinit_constraints();
	system_dof_map.prepare_send_list();
	system_dof_map.compute_sparsity(system_mesh);
};
}

#endif /* COMMON_COMMON_FUNCTIONS_H_ */
