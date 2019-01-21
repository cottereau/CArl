#include "assemble_functions_elasticity_3D.h"
#include "assemble_functions_elasticity_3D_dyn.h"
#include "common_assemble_functions_elasticity_3D.h"


void set_homogeneous_physical_properties_dyn(libMesh::EquationSystems& es, 
        std::string& physicalParamsFile,
        system_weight,system_type,boundary_id_cube,traction_density); 
{
    const libMesh::Parallel::Communicator& SysComm = es.comm();
    int rank = SysComm.rank();
    int nodes = SysComm.size();
     
    WeightFunctionSystemType system_type;
    std::vector<double> force_vol;
    std::vector<double> traction_density;
    int traction_boundary_id;
    // Read the random data info
    double inputE;
    double inputMu;
    
    if(rank == 0)
    {
        std::ifstream physicalParamsIFS(physicalParamsFile);
        physicalParamsIFS >> inputE >> inputMu;
    }

    if(nodes > 1)
    {
        SysComm.broadcast(inputE);
        SysComm.broadcast(inputMu);
    }
    
    // Mesh pointer
    const libMesh::MeshBase& mesh = es.get_mesh();

    // Physical system and its "variables"
    libMesh::ExplicitSystem& physical_param_system = es.get_system<libMesh::ExplicitSystem>("PhysicalConstants");
    const libMesh::DofMap& physical_dof_map = physical_param_system.get_dof_map();
    
    unsigned int physical_consts[2];
    physical_consts[0] = physical_param_system.variable_number ("E");
    physical_consts[1] = physical_param_system.variable_number ("mu");
    
    std::vector<libMesh::dof_id_type> physical_dof_indices_var;
    libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    
    for ( ; el != end_el; ++el)
    {
        const libMesh::Elem* elem = *el;
        // Young modulus, E
        physical_dof_map.dof_indices(elem, physical_dof_indices_var, physical_consts[0]);
        libMesh::dof_id_type dof_index = physical_dof_indices_var[0];
        
        if( (physical_param_system.solution->first_local_index() <= dof_index) &&
        (dof_index < physical_param_system.solution->last_local_index()) )
        {
            physical_param_system.solution->set(dof_index, inputE);
        }
   
        // Mu
        physical_dof_map.dof_indices (elem, physical_dof_indices_var, physical_consts[1]);
        
        dof_index = physical_dof_indices_var[0];
        
        if( (physical_param_system.solution->first_local_index() <= dof_index) &&
        (dof_index < physical_param_system.solution->last_local_index()) )
        {
            physical_param_system.solution->set(dof_index, inputMu);
        }
    }
    physical_param_system.solution->close();
    physical_param_system.update();
}



void set_heterogeneous_physical_properties_dyn(libMesh::EquationSystems& es, std::string& physicalParamsFile)
{
	const libMesh::Parallel::Communicator& SysComm = es.comm();
	int rank = SysComm.rank();
	int nodes = SysComm.size();

	// Read the random data info
	std::vector<double> inputE;
	std::vector<double> inputMu;
	std::vector<double> inputRho;
	std::vector<int> 	inputIdx;

	int NbOfSubdomains = -1;

	double meanE = 0;
	double meanMu = 0;
	double meanRho = 0;

	if(rank == 0)
	{
		std::ifstream physicalParamsIFS(physicalParamsFile);
		physicalParamsIFS >> NbOfSubdomains;
		inputE.resize(NbOfSubdomains);
		inputMu.resize(NbOfSubdomains);
		inputRho.resize(NbOfSubdomains);
		inputIdx.resize(NbOfSubdomains);

		for(int iii = 0; iii < NbOfSubdomains; ++iii)
		{
			physicalParamsIFS >> inputE[iii];
			physicalParamsIFS >> inputMu[iii];
			physicalParamsIFS >> inputRho[iii];
			physicalParamsIFS >> inputIdx[iii];

			meanE += inputE[iii];
			meanMu += inputMu[iii];
			meanRho += inputRho[iii];
		}
		meanE /= NbOfSubdomains;
		meanMu /= NbOfSubdomains;
		meanRho /= NbOfSubdomains;
		physicalParamsIFS.close();
	}

	if(nodes > 1)
	{
		SysComm.broadcast(NbOfSubdomains);
		SysComm.broadcast(meanE);
		SysComm.broadcast(meanMu);
		SysComm.broadcast(meanRho);

		if(rank != 0)
		{
			inputE.resize(NbOfSubdomains);
			inputMu.resize(NbOfSubdomains);
			inputRho.resize(NbOfSubdomains);
			inputIdx.resize(NbOfSubdomains);
		}
		SysComm.broadcast(inputE);
		SysComm.broadcast(inputMu);
		SysComm.broadcast(inputRho);
		SysComm.broadcast(inputIdx);
	}

	// Mesh pointer
	const libMesh::MeshBase& mesh = es.get_mesh();

	// Physical system and its "variables"
	libMesh::ExplicitSystem& physical_param_system = es.get_system<libMesh::ExplicitSystem>("PhysicalConstants");
	const libMesh::DofMap& physical_dof_map = physical_param_system.get_dof_map();

	unsigned int physical_consts[3];
	physical_consts[0] = physical_param_system.variable_number ("E");
	physical_consts[1] = physical_param_system.variable_number ("mu");
	physical_consts[2] = physical_param_system.variable_number ("rho");

	std::vector<libMesh::dof_id_type> physical_dof_indices_var;
	libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	int currentSubdomain = -1;

	for ( ; el != end_el; ++el)
	{
		const libMesh::Elem* elem = *el;

		currentSubdomain = elem->subdomain_id()-1;

		// Young modulus, E
		physical_dof_map.dof_indices(elem, physical_dof_indices_var, physical_consts[0]);
		libMesh::dof_id_type dof_index = physical_dof_indices_var[0];

		if( (physical_param_system.solution->first_local_index() <= dof_index) &&
		(dof_index < physical_param_system.solution->last_local_index()) )
		{
			physical_param_system.solution->set(dof_index, inputE[currentSubdomain]);
		}

		// Mu
		physical_dof_map.dof_indices (elem, physical_dof_indices_var, physical_consts[1]);

		dof_index = physical_dof_indices_var[0];

		if( (physical_param_system.solution->first_local_index() <= dof_index) &&
		(dof_index < physical_param_system.solution->last_local_index()) )
		{
			physical_param_system.solution->set(dof_index, inputMu[currentSubdomain]);
		}
		// Rho
		physical_dof_map.dof_indices (elem, physical_dof_indices_var, physical_consts[2]);

		dof_index = physical_dof_indices_var[1];

		if( (physical_param_system.solution->first_local_index() <= dof_index) &&
		(dof_index < physical_param_system.solution->last_local_index()) )
		{
			physical_param_system.solution->set(dof_index, inputRho[currentSubdomain]);
		}
	}

	physical_param_system.solution->close();
	physical_param_system.update();
}


/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=4 et tw=80 smartindent :                               */

