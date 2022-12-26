/*
 * newmark_param_parser.cpp
 *
 *  Created on: August 3, 2022
 *      Author: Chensheng Luo
 */

#include "newmark_param_parser.h"

namespace carl
{

void get_newmark_params(GetPot& field_parser,
    NewmarkParams& newmark) {

  std::cout << "-- Begin reading Newmark parameters......" << std::endl;

  if (field_parser.search(3, "deltat","Deltat","DELTAT")) {
    newmark.deltat = field_parser.next(newmark.deltat);
    std::cout << "-- Delta t: "<< newmark.deltat << std::endl;
  }else{
    newmark.deltat = 0.001;
    std::cout << "-- Delta t: "<< newmark.deltat << "[No entry, default value is taken]" << std::endl;
  }
  

  if (field_parser.search(3, "alpha","Alpha","ALPHA")) {
    newmark.alpha = field_parser.next(newmark.alpha);
    std::cout << "-- Alpha: "<< newmark.alpha << std::endl;
  }else{
    newmark.alpha = 0;
    std::cout << "-- Alpha: "<< newmark.alpha << "[No entry, default value is taken]"<< std::endl;
  }
  

  if (field_parser.search(3, "beta","Beta","BETA")) {
    newmark.beta = field_parser.next(newmark.beta);
    std::cout << "-- Beta: "<< newmark.beta << std::endl;
  }else{
    newmark.beta = 0.25;
    std::cout << "-- Beta: "<< newmark.beta << "[No entry, default value is taken]"<< std::endl;
  }
  
  if (field_parser.search(3, "gamma","Gamma","GAMMA")) {
    newmark.gamma = field_parser.next(newmark.gamma);
    std::cout << "-- Gamma: "<< newmark.gamma << std::endl;
  }else{
    newmark.gamma = 0.5;
    std::cout << "-- Gamma: "<< newmark.gamma << "[No entry, default value is taken]"<< std::endl;
  }

  std::cout << "...... Reading Newmark parameters FINISH!" << std::endl;
  
};
}