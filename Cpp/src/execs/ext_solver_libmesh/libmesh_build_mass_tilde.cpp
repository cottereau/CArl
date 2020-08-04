// This example generate coupling matrix
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>

#include "libmesh_assemble_lin_homogeneous.h"

using namespace libMesh;

template <class Container> void split(const std::string& str, Container& cont, char delim = ' ')
{
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        cont.push_back(token);
    }
}


int main(int argc, char** argv) {
	 // --- Initialize libMesh
libMesh::LibMeshInit init(argc, argv);

  // libMesh's C++ / MPI communicator wrapper
libMesh::Parallel::Communicator& WorldComm = init.comm();

  // Number of processors and processor rank.
int rank = WorldComm.rank();
int nodes = WorldComm.size();

  // Command line parser
GetPot command_line(argc, argv);

  // File parser
GetPot field_parser;

  // If there is an input file, parse it to get the parameters. Else, parse the command line
std::string input_filename;
if (command_line.search(2, "--inputfile", "-i")) {
	input_filename = command_line.next(input_filename);
    field_parser.parse_input_file(input_filename, "#", "\n", " \t\n");
} else{
    field_parser = command_line;
}

  //Parse file content to the input_params 
libmesh_assemble_input_params input_params;
get_input_params(field_parser, input_params);

  //Preparation to build the PETSC matrix M, K
Mat sys_mat_PETSC_M;
Mat sys_mat_PETSC_K;

libMesh::Real delta_t = input_params.deltat;
libMesh::Real beta = 0.25;
PetscScalar a = beta*delta_t*delta_t;

MatCreate(WorldComm.get(),&sys_mat_PETSC_M);
MatCreate(WorldComm.get(),&sys_mat_PETSC_K);

carl::read_PETSC_matrix(sys_mat_PETSC_M, input_params.output_base+"_mass_test_sys_mat.petscmat", WorldComm.get());
carl::read_PETSC_matrix(sys_mat_PETSC_K, input_params.output_base+"_stiffness_test_sys_mat.petscmat", WorldComm.get());

//calculate de M~  = M + a*K
MatAXPY(sys_mat_PETSC_M, a, sys_mat_PETSC_K, DIFFERENT_NONZERO_PATTERN);

libMesh::PetscMatrix<libMesh::Number> sys_mat_PETSC_M_tilde(sys_mat_PETSC_M,WorldComm);


std::string mass_tilde_path = "./GC_solver/mass_tilde/"; 
std::vector<std::string> str;
split(input_params.output_base, str,'/');

std::string name_matrix_tilde_matlab = mass_tilde_path+"mass_tilde_"+str[str.size()-1]+".m";
std::string name_matrix_tilde_petsc  = mass_tilde_path+"mass_tilde_"+str[str.size()-1]+".petscmat";

#ifdef PRINT_MATLAB_DEBUG
  sys_mat_PETSC_M_tilde.print_matlab(name_matrix_tilde_matlab);
#endif

carl::write_PETSC_matrix(sys_mat_PETSC_M_tilde,name_matrix_tilde_petsc);


	return 0;
}


