/*
 * weight_parameter_function.h
 *
 *  Created on: Feb 2, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef WEIGHT_PARAMETER_FUNCTION_H_
#define WEIGHT_PARAMETER_FUNCTION_H_

#include "common_header_ext_solver_libmesh.h"

class weight_parameter_function
{
protected:
  // Members
//  const libMesh::PointLocatorBase& m_locator;
  libMesh::Mesh& m_alpha_mesh;
  std::unique_ptr<libMesh::PointLocatorBase> m_locator_unique_ptr;

  double m_alpha_eps;
  double m_alpha_coupling_BIG;
  double m_alpha_coupling_micro;

  int m_subdomain_idx_BIG;
  int m_subdomain_idx_micro;
  int m_subdomain_idx_both;
  int m_subdomain_idx_coupling;

public:

  // Constructor
  weight_parameter_function(  libMesh::Mesh& alpha_mesh, double alpha_eps,
            double alpha_coupling_BIG,
            int subdomain_idx_BIG, int subdomain_idx_micro, int subdomain_idx_both, int subdomain_idx_coupling) :
    m_alpha_mesh{ alpha_mesh },
    m_locator_unique_ptr { alpha_mesh.sub_point_locator() },
    m_alpha_eps { alpha_eps },
    m_alpha_coupling_BIG { alpha_coupling_BIG },
    m_alpha_coupling_micro { 1. - alpha_coupling_BIG },
    m_subdomain_idx_BIG { subdomain_idx_BIG },
    m_subdomain_idx_micro { subdomain_idx_micro },
    m_subdomain_idx_both { subdomain_idx_both },
    m_subdomain_idx_coupling { subdomain_idx_coupling }
  {
  };

  weight_parameter_function(  libMesh::Mesh& alpha_mesh ) :
    m_alpha_mesh{ alpha_mesh },
    m_locator_unique_ptr { alpha_mesh.sub_point_locator() },
    m_alpha_eps { 10E-2 },
    m_alpha_coupling_BIG { 0.5 },
    m_alpha_coupling_micro { 0.5 },
    m_subdomain_idx_BIG { -1 },
    m_subdomain_idx_micro { -1 },
    m_subdomain_idx_both { -1 },
    m_subdomain_idx_coupling { -1 }
  {
  };

  // Destructor
  ~weight_parameter_function()
  {
    clear();
  }

  void set_parameters(double alpha_eps, double alpha_coupling_BIG,
            int subdomain_idx_BIG, int subdomain_idx_micro, int subdomain_idx_both, int subdomain_idx_coupling)
  {
    m_alpha_eps = alpha_eps;

    m_alpha_coupling_BIG = alpha_coupling_BIG;
    m_alpha_coupling_micro = 1. - alpha_coupling_BIG;

    m_subdomain_idx_BIG  = subdomain_idx_BIG;
    m_subdomain_idx_micro = subdomain_idx_micro;
    m_subdomain_idx_both = subdomain_idx_both;
    m_subdomain_idx_coupling = subdomain_idx_coupling;
  };

  void set_parameters( std::string &filename )
  {
    GetPot field_parser;
    field_parser.parse_input_file(filename, "#", "\n", " \t\n");

    if (field_parser.search(1, "MacroDomainOnlyIdx")) {
      m_subdomain_idx_BIG = field_parser.next(
          m_subdomain_idx_BIG);
    } else {
      homemade_error_msg("Missing the macro domain only index!");
    }

    if (field_parser.search(1, "MicroDomainOnlyIdx")) {
      m_subdomain_idx_micro = field_parser.next(
          m_subdomain_idx_micro);
    } else {
      homemade_error_msg("Missing the micro domain only index!");
    }

    if (field_parser.search(1, "CouplingDomainIdx")) {
      m_subdomain_idx_coupling = field_parser.next(
          m_subdomain_idx_coupling);
    } else {
      homemade_error_msg("Missing the coupling domain index!");
    }

    if (field_parser.search(1, "BothDomainsIdx")) {
      m_subdomain_idx_both = field_parser.next(
          m_subdomain_idx_both);
    } else {
      m_subdomain_idx_both = -1;
      std::cout << " >> Warning: Both system's domain index not defined!" << std::endl;
    }

    if (field_parser.search(1, "AlphaCouplingMacro")) {
      m_alpha_coupling_BIG = field_parser.next(
          m_alpha_coupling_BIG);
      m_alpha_coupling_micro = 1. - m_alpha_coupling_BIG;
    } else if (field_parser.search(1, "AlphaCouplingMicro")) {
      m_alpha_coupling_micro = field_parser.next(
          m_alpha_coupling_micro);
      m_alpha_coupling_BIG = 1. - m_alpha_coupling_micro;
    } else {
      m_alpha_coupling_BIG = 0.5;
      m_alpha_coupling_micro = 0.5;
      std::cout << " >> Warning: Coupling weights not defined, using 0.5 for both systems!" << std::endl;
    }

    if (field_parser.search(1, "AlphaEps")) {
      m_alpha_eps = field_parser.next(m_alpha_eps);
    } else {
      m_alpha_eps = 1e-2;
      std::cout << " >> Warning: Alpha epsilon not defined, using 1e-2!" << std::endl;
    }
  }

  double get_alpha(const libMesh::Point& qpoint, WeightFunctionSystemType system_type)
  {
    switch (system_type)
    {
      case WeightFunctionSystemType::MICRO :  
          return this->get_alpha_micro(qpoint);
          break;
      case WeightFunctionSystemType::MACRO :  
          return this->get_alpha_BIG(qpoint);
          break;
      case WeightFunctionSystemType::NO_WEIGHT :  
          return 1.;
          break;
    }
    homemade_error_msg("Why are you here!?");
  }

  double get_alpha_BIG(const libMesh::Point& qpoint)
  {
    libMesh::PointLocatorBase& locator = *m_locator_unique_ptr.get();
    const libMesh::Elem* elem = locator(qpoint);
    double output = 0;

    if(elem->subdomain_id() == m_subdomain_idx_BIG)
    {
      output = 1;
    }
    else if(elem->subdomain_id() == m_subdomain_idx_micro)
    {
      output = 0;
    }
    else if(elem->subdomain_id() == m_subdomain_idx_both)
    {
      output = m_alpha_eps;
    }
    else if(elem->subdomain_id() == m_subdomain_idx_coupling)
    {
      output = m_alpha_coupling_BIG;
    }

    return output;
  }

  double get_alpha_micro(const libMesh::Point& qpoint)
  {
    libMesh::PointLocatorBase& locator = *m_locator_unique_ptr.get();
    const libMesh::Elem* elem = locator(qpoint);
    double output = 0;

    if(elem->subdomain_id() == m_subdomain_idx_BIG)
    {
      output = 0;
    }
    else if(elem->subdomain_id() == m_subdomain_idx_micro)
    {
      output = 1;
    }
    else if(elem->subdomain_id() == m_subdomain_idx_both)
    {
      output = 1 - m_alpha_eps;
    }
    else if(elem->subdomain_id() == m_subdomain_idx_coupling)
    {
      output = m_alpha_coupling_micro;
    }

    return output;
  }

  // Methods
  void clear()
  {
    m_locator_unique_ptr.reset(NULL);
  };
};

#endif /* WEIGHT_PARAMETER_FUNCTION_H_ */
