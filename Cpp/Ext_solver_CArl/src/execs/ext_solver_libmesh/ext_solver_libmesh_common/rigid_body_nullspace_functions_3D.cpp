#include "rigid_body_nullspace_functions_3D.h"

void build_rigid_body_vectors(libMesh::ImplicitSystem&  input_system, MatNullSpace& nullsp_sys)
{
	// --- Set up some temporary variables to simplify code

	// System matrix pointer (only used to get the proper vector dimensions)
	libMesh::PetscMatrix<libMesh::Number> * mat_sys = libMesh::cast_ptr<libMesh::PetscMatrix<libMesh::Number>* >(input_system.matrix);

	/// Mesh address
	const libMesh::MeshBase& mesh_sys = input_system.get_mesh();

	// System number
	unsigned int sys_number = input_system.number();

	// --- Set the coordinates vector structure as the same of rigidity matrix
	Vec coord_vec_PETSC;
	PetscInt local_N;
	MatGetLocalSize(mat_sys->mat(),NULL,&local_N);
	VecCreate(mesh_sys.comm().get(),&coord_vec_PETSC);
	VecSetSizes(coord_vec_PETSC,local_N,mat_sys->n());
	VecSetBlockSize(coord_vec_PETSC,mesh_sys.mesh_dimension());
	VecSetFromOptions(coord_vec_PETSC);

	libMesh::PetscVector<libMesh::Number> coord_vec(coord_vec_PETSC,mesh_sys.comm());

	// --- Fill the coordinates vector
	auto node_it = mesh_sys.local_nodes_begin();
	auto node_it_end = mesh_sys.local_nodes_end();

	unsigned int dof_number = 0;

	for( ; node_it != node_it_end; ++node_it)
	{
		const libMesh::Node* node = *node_it;

		for(unsigned int var=0; var<node->n_dofs(sys_number); var++)
		{
			dof_number = node->dof_number(sys_number,var,0);
			coord_vec.set(dof_number,node->operator ()(var));
		}
	}

	MatNullSpaceCreateRigidBody(coord_vec.vec(),&nullsp_sys);
	VecDestroy(&coord_vec_PETSC);
};


void write_rigid_body_vectors(MatNullSpace& nullsp_sys, const std::string output_base, int rank)
{
	// PETSc variables
	PetscInt    nullsp_nvecs;
	const Vec * nullsp_vecs;
	PetscBool   nullsp_has_const;
	
	// Get the vectors
	MatNullSpaceGetVecs(nullsp_sys,&nullsp_has_const,&nullsp_nvecs,&nullsp_vecs);

	// Print them
	std::string filename;
	for(int iii = 0; iii < nullsp_nvecs; ++iii)
	{
		filename = output_base + "_rb_vector_" + std::to_string(iii) + "_n_" + std::to_string(nullsp_nvecs);
		carl::write_PETSC_vector_MATLAB(nullsp_vecs[iii],filename + ".m",PETSC_COMM_WORLD);

		// Export the RB vector
		carl::write_PETSC_vector(nullsp_vecs[iii], filename + ".petscvec",rank,PETSC_COMM_WORLD);
	};
};
