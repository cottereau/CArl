// This example generate interpolation matrix
#include "CArl_build_interpolation.h"
/** \brief Program used to calculate interpolation matrix.
 * 
 *  Usage: `./libmesh_build_interpolation_matrix -i [input file]`
 *  
 * The input file is parsed by the get_input_params(GetPot& field_parser, libmesh_assemble_input_params& input_params) function, and it contains the following parameters. 
 *
 *  Required parameters:
 *    - 'AlphaA' : gamma_A*delta_T, the product of the Newmark coefficient and the time interval for part A
 *    - 'AlphaB' : gamma_B*delta_t, the product of the Newmark coefficient and the time interval for part B
 *    - 'mTildeA' : Sum of Mass matrix and Newmark matrix for part A
 *    - 'mTildeB' : Sum of Mass matrix and Newmark matrix for part B
 *    - 'couplingMatrixA' : Coupling matrix for part A
 *    - 'couplingMatrixB' : Coupling matrix for part B
 *    - 'OutputFolder' : Folder where the result is stored
 *    - 'InverseNameA' : Name in which the inverse matrix of Mass for part A is stored 
 *    - 'InverseNameB' : Name in which the inverse matrix of Mass for part B is stored 
 */
using namespace std;
using namespace libMesh;

void load_PETSC_matrices(std::string path_dynamic_matrices, 
        std::string path_coupling_matix,
        libMesh::Parallel::Communicator& WorldComm,
        Mat& M_tilde,
        Mat& dense_M_tilde,
        Mat& Coupling);

void new_invert_PETSC_matrix(libMesh::Parallel::Communicator& WorldComm,
        Mat& dense_matrix,
        Mat &inv_matrix,
        std::string inv_path);

int main(int argc, char** argv) 
{
	 // --- Initialize libMesh
libMesh::LibMeshInit init(argc, argv);

  // libMesh's C++ / MPI communicator wrapper
libMesh::Parallel::Communicator& WorldComm = init.comm();

  // Number of processors and processor rank.
int rank = WorldComm.rank();
int nodes = WorldComm.size();

  // --- Set up input_params
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

// Parse file content to the input_params
carl_build_interpolation_params input_params;
get_input_params(field_parser, input_params);

// --- Definition of matrix
// Interpolation matrix
Mat* sys_mat_PETSC_H_save;

// Mass matrix and its dense matrix(the inverse of mass matrix is done by the dense matrix)
Mat m_tilde_A;
Mat m_tilde_B;
Mat dense_m_tilde_A;
Mat dense_m_tilde_B;

// Inversed matrix
Mat inv_m_tilde_A;
Mat inv_m_tilde_B;

Mat coupling_matrix_A;
Mat coupling_matrix_B;

Mat H_1;
Mat H_2;
PetscErrorCode ierr;
PetscScalar alpha_A = input_params.alpha_A;
PetscScalar alpha_B = input_params.alpha_B;

std::cout << " before load_PETSC_matrices " << std::endl;//debug usage

load_PETSC_matrices(input_params.path_tilde_matrix_A,
              input_params.path_macro_coupling_matrix,
              WorldComm,
              m_tilde_A, 
              dense_m_tilde_A,
              coupling_matrix_A);

load_PETSC_matrices(input_params.path_tilde_matrix_B,
              input_params.path_micro_coupling_matrix,
              WorldComm,
              m_tilde_B,
              dense_m_tilde_B,
              coupling_matrix_B);

std::cout << " after load_PETSC_matrices " << std::endl;//debug usage

// --- Compute part

// - Invert mass matrix
new_invert_PETSC_matrix(WorldComm,
              dense_m_tilde_A, 
              inv_m_tilde_A,
              input_params.output_base+input_params.path_inv_matrix_A);

new_invert_PETSC_matrix(WorldComm,
              dense_m_tilde_B, 
              inv_m_tilde_B,
              input_params.output_base+input_params.path_inv_matrix_B);

// - Product coupling and mass matrix

// Examination of dimension
PetscInt        petsc_m=0, petsc_n=0,petsc_a=0,petsc_b=0,petsc_c=0,petsc_d=0;
MatGetSize(coupling_matrix_A, &petsc_m, &petsc_n);
MatGetSize(inv_m_tilde_A, &petsc_a, &petsc_b);
//MatGetSize(t_coupling_matrix_A, &petsc_c, &petsc_d);
if (petsc_a != petsc_n || petsc_b != petsc_n){
  homemade_error_msg("Dimension Error! The coupling matrix and the mass matrix are not compatible for matrix product");
}

// Create result matrix
MatCreateAIJ(WorldComm.get(),PETSC_DECIDE,PETSC_DECIDE,petsc_m,petsc_m,petsc_m,NULL,petsc_m,NULL,&H_1);
MatCreateAIJ(WorldComm.get(),PETSC_DECIDE,PETSC_DECIDE,petsc_m,petsc_m,petsc_m,NULL,petsc_m,NULL,&H_2);


//H_1=C_A*M_A*C_A^T
MatRARt(inv_m_tilde_A, 
        coupling_matrix_A,
        MAT_INITIAL_MATRIX,
        PETSC_DEFAULT,
        &H_1);

//H_2=C_B*M_B*C_B^T
MatRARt(inv_m_tilde_B, 
        coupling_matrix_B,
        MAT_INITIAL_MATRIX,
        PETSC_DEFAULT,
        &H_2);
              

// - Linear scale and  addition
//H_2=H_2*alpha_B
MatScale(H_2,alpha_B);

//H_1=H_1*alpha_A+H_2
MatAYPX(H_1,alpha_A,H_2, SAME_NONZERO_PATTERN);// Not sure on the MatStructure str
sys_mat_PETSC_H_save = &H_1;


std::cout << " Calculate finish " << std::endl;//debug usage


// --- Output interpolation matrix
libMesh::PetscMatrix<libMesh::Number> sys_H_mat(*sys_mat_PETSC_H_save,WorldComm);
std::cout << " Output Interpolation Matrix " << std::endl;//debug usage

// Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
   sys_H_mat.print_matlab(input_params.output_base + "_interpolation_mat.m");
#endif

carl::write_PETSC_matrix(*sys_mat_PETSC_H_save, input_params.output_base + "_interpolation_mat.petscmat",0,WorldComm.get(),1);

return 0;
}

void load_PETSC_matrices(std::string path_dynamic_matrices, 
        std::string path_coupling_matix,
        libMesh::Parallel::Communicator& WorldComm,
        Mat& M_tilde,
        Mat& dense_M_tilde,
        Mat& Coupling)
{
  // Create matrix
  MatCreate(PETSC_COMM_SELF,&dense_M_tilde);
  MatSetType(dense_M_tilde,MATSEQAIJ);
  MatCreate(WorldComm.get(),&M_tilde);
  MatCreate(WorldComm.get(),&Coupling);


  // Read Matrix
  carl::read_PETSC_matrix(Coupling, 
                path_coupling_matix, WorldComm.get());
  carl::read_PETSC_matrix(M_tilde, 
                path_dynamic_matrices, WorldComm.get());
  carl::read_PETSC_matrix(dense_M_tilde, 
                path_dynamic_matrices,PETSC_COMM_SELF);
  
}

void new_invert_PETSC_matrix(libMesh::Parallel::Communicator& WorldComm,
          Mat& dense_matrix, 
          Mat& inv_matrix,
          std::string inv_path)
{
  Mat             Ones,inv_dense_matrix;
  PetscInt        i,j;
  PetscErrorCode  ierr=0;
  PetscScalar     v = 1.0;
  MatFactorInfo   factor_info;
  IS        rperm, cperm;
  const PetscInt    *cols;
  const PetscScalar *vals;

  //Create identity and inverse of the dense matrix
  MatDuplicate(dense_matrix,MAT_DO_NOT_COPY_VALUES,&inv_dense_matrix);
  MatDuplicate(dense_matrix,MAT_DO_NOT_COPY_VALUES,&Ones);
  MatSetType(inv_dense_matrix, MATSEQDENSE);
  MatSeqDenseSetPreallocation(inv_dense_matrix,NULL);
  MatSetType(Ones, MATSEQDENSE);
  MatSeqDenseSetPreallocation(Ones,NULL);
  MatZeroEntries(Ones);
  MatShift(Ones,1);
  
  // Factor input matrix
  MatGetOrdering(dense_matrix,MATORDERINGNATURAL,  &rperm,  &cperm);
  MatFactorInfoInitialize(&factor_info);
  MatLUFactor(dense_matrix,rperm,cperm,&factor_info);
  
  // Calculate inverse of the dense matrix 
  ierr = MatMatSolve(dense_matrix,Ones,inv_dense_matrix);

  // Reset input's factoring
  MatSetUnfactored(dense_matrix);

  //Output inverse of the dense matrix
  carl::write_PETSC_matrix(inv_dense_matrix,inv_path+".petscmat",0,PETSC_COMM_SELF,1);

  #ifdef PRINT_MATLAB_DEBUG
    carl::write_PETSC_matrix_MATLAB(inv_dense_matrix,inv_path+".m",PETSC_COMM_SELF);
  #endif

  // Cleanup
  MatDestroy(&Ones);
  MatDestroy(&dense_matrix);
  MatDestroy(&inv_dense_matrix);

  //Load Matrix in WorldComm
  MatCreate(WorldComm.get(),&inv_matrix);//changed--> changement!
  carl::read_PETSC_matrix(inv_matrix, 
                inv_path+".petscmat", WorldComm.get());
}

