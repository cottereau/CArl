// This example generate interpolation matrix
#include "libmesh_assemble_lin_homogeneous.h"

using namespace libMesh;

void load_PETSC_matrices(std::string path_dynamic_matrices, 
        std::string path_coupling_matix,
        libMesh::Parallel::Communicator& WorldComm,
        Mat& M_tilde,
        Mat& Coupling);

void inv_PETSC_matrix(Mat& matrix, 
        Mat &inv_matrix);

int main(int argc, char** argv) 
{
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
if (command_line.search(2, "--inputfile", "-i")) 
{
  	input_filename = command_line.next(input_filename);
    field_parser.parse_input_file(input_filename, "#", "\n", " \t\n");
} 
else
{
  field_parser = command_line;
}

//Parse file content to the input_params
libmesh_assemble_interpolation_params input_params;
//	get_input_params(field_parser, input_params);

//Interpolation matrix
Mat sys_mat_PETSC_H_save;

// Load all the matrices of A domain needed
Mat m_tilde_A;
Mat coupling_matrix_A;

// Transformed matrix
Mat inv_m_tilde_A;
Mat t_coupling_matrix_A;

Mat* H_1;
libMesh::Real alpha_A = input_params.alpha_A;

printf("befor load_PETSC_matrices");

load_PETSC_matrices(input_params.path_tilde_matrices_A,
              input_params.path_macro_coupling_matrix,
              WorldComm,
              m_tilde_A, 
              coupling_matrix_A);

printf("after load_PETSC_matrices\n");

//Compute part

inv_PETSC_matrix(m_tilde_A, inv_m_tilde_A);

MatTranspose(coupling_matrix_A, MAT_INPLACE_MATRIX, &t_coupling_matrix_A);

MatMatMatMult(coupling_matrix_A, 
              inv_m_tilde_A, 
              t_coupling_matrix_A, 
              MAT_REUSE_MATRIX,  
              PETSC_DEFAULT, 
              H_1);

// // Print MatLab debugging output? Variable defined at "carl_headers.h"
// #ifdef PRINT_MATLAB_DEBUG
//   sys_mat_PETSC_H_save.print_matlab(input_params.output_folder + "m_tilde_A.m");
// #endif

	return 0;
}

void load_PETSC_matrices(std::string path_dynamic_matrices, 
        std::string path_coupling_matix,
        libMesh::Parallel::Communicator& WorldComm,
        Mat& M_tilde,
        Mat& Coupling)
{

  MatCreate(WorldComm.get(),&M_tilde);
  MatCreate(WorldComm.get(),&Coupling);

  carl::read_PETSC_matrix(M_tilde, 
                path_dynamic_matrices, WorldComm.get());

  carl::read_PETSC_matrix(Coupling, 
                path_coupling_matix, WorldComm.get());
}

void inv_PETSC_matrix(Mat& matrix,  Mat &inv_matrix)
{

  Mat             Ones;
  PetscInt        petsc_m=0, petsc_n=0,i,j;
  PetscErrorCode  ierr=0;
   PetscScalar    v = 1.0;  
  ierr = MatGetSize(matrix, &petsc_m, &petsc_n);

    //Create a Ones matrix
  MatCreateSeqDense(PETSC_COMM_SELF,petsc_m,petsc_n,NULL,&Ones);
    //Create a Solution matrix
  MatCreateSeqDense(PETSC_COMM_SELF,petsc_m,petsc_n,NULL,&inv_matrix);

  //Factorisation LU
  MatLUFactor(matrix,NULL,NULL,NULL);

    //Affect value to the matrix Ones
  for(i = 0; i<petsc_m; i++)
    for(j=0; j<petsc_n; j++)
      MatSetValues(Ones,1,&i,1,&j,&v,INSERT_VALUES);
  //Inversion
  MatMatSolve(matrix, Ones, inv_matrix);
}