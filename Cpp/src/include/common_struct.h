/*
 * common_struct.h
 *
 *  Created on: June 25, 2022
 *      Author: Chensheng Luo
 */

#ifndef COMMON_STRUCT_H_
#define COMMON_STRUCT_H_

namespace carl
{

// [DYN]
struct DynInitialVectorPath {
    std::string disp;
    std::string speed;
};
// [DYN]
struct DynSystemMatrixPath {
    std::string mass_tilde;          ///< mass tilde matrix
    std::string coupling;       ///< coupling matrix
    std::string stiffness;      ///< stiffness matrix
    std::string damping;       ///< damping matrix
};
// [DYN]
struct DynSystemVectorPath
  {
    std::string prev_acc;
    std::string prev_disp;
    std::string prev_disp_free;
    std::string prev_speed;
    std::string prev_speed_free;

    std::string rhs_free;
    std::string rhs_link;

    std::string this_acc;
    std::string this_acc_free;
    std::string this_acc_link;
    std::string this_disp;
    std::string this_disp_free;
    std::string this_disp_link;
    std::string this_speed;
    std::string this_speed_free;
    std::string this_speed_link;

    std::string inter_disp_free;

    std::string this_force;
    std::string next_force;
    int coupling_sign;
  };

}





#endif /* COMMON_STRUCT_H_ */