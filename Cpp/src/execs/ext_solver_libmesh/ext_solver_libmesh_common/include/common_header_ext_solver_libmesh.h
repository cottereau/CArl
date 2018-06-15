/*
 * common_est_libmesh_header.h
 *
 *  Created on: Apr 11, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_EXT_LIBMESH_HEADER_H_
#define COMMON_EXT_LIBMESH_HEADER_H_

#include "libmesh/analytic_function.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/perf_log.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/point_locator_tree.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_tetgen_interface.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/metis_partitioner.h"
#include "libmesh/error_vector.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/statistics.h"
#include "libmesh/namebased_io.h"
#include "libmesh/petsc_preconditioner.h"
#include "libmesh/fem_system.h"

#include "libmesh/boundary_info.h"
#include "libmesh/diff_solver.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_base.h"
#include "libmesh/fem_context.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/enum_solver_type.h"
#include "libmesh/newton_solver.h"
#include "libmesh/newmark_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/steady_solver.h"
#include "libmesh/transient_system.h"

#include "libmesh/dense_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"
#include "libmesh/getpot.h"

#include "libmesh/nonlinear_solver.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/nonlinear_implicit_system.h"

#include <petscmat.h>
#include <petscvec.h>
#include <petscsys.h>
#include <petscksp.h>

#include "carl_headers.h"
// #define homemade_error_msg(msg) do \
// { \
// 	std::cerr << "Error: " <<  msg << std::endl; \
// 	std::exit(EXIT_FAILURE); \
// } while(false)

#endif /* COMMON_EXT_LIBMESH_HEADER_H_ */
