/*
 * carl_loop_dyn_input_parser.h 
 *
 *  Created on: Nov 17, 2021
 *      Author: Chensheng Luo
 */

#ifndef CARL_LOOP_DYN_INPUT_PARSER_H_
#define CARL_LOOP_DYN_INPUT_PARSER_H_

#include "carl_headers.h"
#include "newmark_param_parser.h"
#include <algorithm>

namespace carl
{
/// Structure containing the parameters for begining the dynamic loop

struct feti_loop_dyn_params {
	// --- Parameters used directly by the CArl_loop_dyn program

	// Cluster 
	ClusterSchedulerType scheduler;     ///< Cluster scheduler software type. *Values*: PBS, SLURM (code not implemented for the later yet).
    std::string script_filename;         ///< Cluster script file name
	
	// Path to "scratch" folder
	std::string scratch_folder_path;	///< Path to the folder which will be used to save the temporary files during the solve operation
    std::string result_folder_path;     ///< Path to the folder which will be used to save the result of all calculation.
	
    // Path of calculated matrix
    struct carl::DynSystemMatrixPath matrix_A; ///< All matrices' path of object A
    struct carl::DynSystemMatrixPath matrix_B; ///< All matrices' path of object B
    std::string coupling_folder_path;    ///< Coupling Matrix folder path
    std::string interpolation_matrix;    ///< Interpolation matrix

    // Initial files
    struct carl::DynInitialVectorPath initial_A; ///< All initial conditions' path of object A
    struct carl::DynInitialVectorPath initial_B; ///< All initial conditions' path of object B

    // Force files
    int force_prepare_method;                    ///< Force preparation method
    std::string force_prepare_params;            ///< Force preparation input file

    // External Solver files
    std::string ext_solver_A_input;                ///< Input parameters of the external solver for system A.
    std::string ext_solver_B_input;                ///< Input parameters of the external solver for system B.
    std::string ext_solver_general_input;          ///< Input parameters of the external solver for general system(to solve interpolation).

    std::string ext_solver_launch_script_A;       ///< Input parameters of the external solve for system A.
    std::string ext_solver_launch_script_B;       ///< Input parameters of the external solve for system A.
    std::string ext_solver_launch_script_general;   ///< Input parameters of the external solve for system A.
    std::string general_entry_file_path;            ///< Path of entry file
    
    // Newmark Parameter
    struct carl::NewmarkParams newmark_A;          ///< All newmark parameters of object A
    struct carl::NewmarkParams newmark_B;          ///< All newmark parameters of object B
    double simulation_duration;                    ///< Simulation duration
    int inner_loop_times;                          ///< Inner loop times
    int outer_loop_times;                          ///< Outer loop times
    int result_times_A;                            ///< Nb of time to print A result
    int result_times_B;                            ///< Nb of time to print B result

    // Rigid body mode options for the micro system
    bool bUseRigidBodyModes;            ///< [RB] Use the rigid body modes for the micro system?
    std::string RB_vectors_base;        ///< [RB] Common path base for the micro system's rigid body mode vectors.
    int nb_of_rb_vectors;               ///< [RB] Number of RB mode vectors.

    int dyn_solver;
    // CG parameters
    double CG_coupled_conv_abs;     ///< [CG] Absolute residual convergence.
    double CG_coupled_conv_rel;     ///< [CG] Relative residual convergence.
    double CG_coupled_div;          ///< [CG] Residual divergence.
    double CG_coupled_conv_corr;    ///< [CG] Relative rigid body mode convergence.
    int CG_coupled_conv_max;        ///< [CG] Maximum number of iterations. 
    carl::BaseCGPrecondType CG_precond_type;    ///< [CG] Type of preconditionner.
};

/** \brief **DYN-DI/DYN-CG** Parser function for dynamic solvers input.
 *  
Parameters:

 + **Scheduler Parameters**: 
     - `ClusterSchedulerType` : the cluster scheduler used in your system, *chosen from* `LOCAL`,`SLURM` or `PBS`   

   *Only applicable for `SLURM` or `PBS`* :   
     - `ScriptFile` : Path to the script file

 + **Solver Parameters**:
     - `DynSolver` : Dynamic solver type,*chosen from* `DI` or `CG`
     - `ScratchFolderPath` : Path to the scratch folder
     - `ResultFolderPath` : Path to the result folder

 + **System Matrix Parameters**:
     - `TildeMatrixA` : Path to the system matrix of A, *i.e.* \f$ \tilde{M}^{A} \f$
     - `StiffnessMatrixA` : Path to the stiffness matrix of A, *i.e.* \f$ K^{A} \f$
     - `DampingMatrixA` : Path to the damping matrix of A, *i.e.* \f$ D^{A} \f$
     - `TildeMatrixB` : Path to the system matrix of B, *i.e.* \f$ \tilde{M}^{B} \f$
     - `StiffnessMatrixB` : Path to the stiffness matrix of B, *i.e.* \f$ K^{B} \f$
     - `DampingMatrixB` : Path to the damping matrix of B, *i.e.* \f$ D^{B} \f$

     - `CouplingMatricesFolder` : Folder to the coupling matrices, generated at the 3rd step(coupling)

     - (**For [CArl-Dyn-DI]**) `InterpolationMatrix` : Path to the interpolation matrix, generated at the 4th step(interpolation) *i.e.* \f$ H \f$
  
  + **Initial Condition Parameters**:
     - (*OPTIONAL*)`InitialDispA` : Path to the initial displacement of A
     - (*OPTIONAL*)`InitialDispB` : Path to the initial displacement of A
     - (*OPTIONAL*)`InitialSpeedA` : Path to the initial speed of A
     - (*OPTIONAL*)`InitialSpeedB` : Path to the initial speed of A     

   * **ATTENTION**: if initial condition is set, do set an equilibrum initial condition!

  + **External Force Parameters**:
     - `ForcePrepareMethod` : Force preparetion method, *chosen from* `ModalSinus`, `ModalConstant`, `ModalLinear`, `ModalProduct`
     - `ForcePrepareParams` : Path to the input file of force preparetion parameters, see get_input_params(GetPot& field_parser,int force_prepare_mode,dyn_force_params& input_params) for its redaction format!
  
  + **External solver parameters**
     - `ExtSolverInputA` : Path to a general external solver input of object A
     - `ExtSolverInputB` : Path to a general external solver input of object B
     - (**For [CArl-Dyn-DI]**)`ExtSolverInputInterpolation` : Path to an external solver input of interpolation matrix/vector   

     * *And either*:   
     - `ExtSolverLaunchScript` : Script to launch external solver for all cases,    

     * *or (if different external solver is needed!)*:   
     - `ExtSolverLaunchScriptA` : Script to launch external solver for object A
     - `ExtSolverLaunchScriptB` : Script to launch external solver for object B
     - (**For [CArl-Dyn-DI]**) `ExtSolverLaunchScriptInterpolation` : Script to launch external solver for interpolation matrix   

   * **Example**:One example of this ExtSolver can be : `'srun -n 4 $CARLBUILD/libmesh_solve_linear_system -i '`
  
  + **Path of this file**:
     - `GeneralEntryFilePath` : Path of this file, which will be served as general extry

  + **Newmark parameters**:      
     - `SimulationDuration` : total duration of simulation, *Default*: 1   
     - `NewmarkParameters` *or* (`NewmarkParametersA` and `NewmarkParametersB`) : Path to the Newmark parameter of two model, see get_newmark_params(GetPot& field_parser,NewmarkParams& newmark) for its redaction format.

   * **Attention**: `SimulationDuration` should be a integer multiple of `deltatA`, whicn should also be a integer multiple of `deltatB`, otherwise error may occur!.
  + (*OPTIONAL WITH DEFAULT*) **Output parameters**:
      In case of very small time step, it's not necessary to show all time step result, so we have following these parameters:

     * *Either*:   
     - `ResultTime` : the time interval to output both result, *Default*: next case

     * *or (if different time interval is needed!)*: 

     - `ResultTimeA` : the time interval to output A's result, *Default*: identical to `deltatA`
     - `ResultTimeB` : the time interval to output B's result, *Default*: identical to `deltatB`

  + (**For [CArl-Dyn-CG]**) **Rigid Body Mode**:
     - Boolean flags `UseRigidBodyModesB` : use the rigid body modes for system B.
     - `ExtForceSystemB` : path to the vector containing the external forces for the system B.
     - `RBVectorBase` : filename base of the rigid body modes vectors.

  + (**For [CArl-Dyn-CG]**, *OPTIONAL WITH DEFAULT*) **CG Solver Parameters**:
     - `CGPreconditionerType` : CG preconditioner type. Values: "NONE", "Coupling_operator" or "Coupling_operator_jacobi". *Default*: "Coupling_operator".
     - `CoupledConvAbs` : CG absolute convergence on the residual. *Default*: 1e-20.
     - `CoupledConvRel` : CG relative convergence on the residual. *Default*: 1e-5.
     - `CoupledCorrConvRel` : CG relative convergence on the rigid body corrections. *Default*: 1e-6.
     - `CoupledDiv` : CG residual divergence parameter. *Default*: 100000.
     - `CoupledIterMax` : CG maximum number of iterations. *Default*: 1000.
 *
 */

void get_input_params(GetPot& field_parser,
		feti_loop_dyn_params& input_params);

bool is_multiple(double big,double small,double tolerance);
};
#endif /* CARL_LOOP_DYN_INPUT_PARSER_H_ */