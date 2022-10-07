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

  if (field_parser.search(3, "deltat","Deltat","DELTAT")) {
    newmark.deltat = field_parser.next(newmark.deltat);
  }
  else{
    newmark.deltat = 0.001;
    std::cout << "[Newmark Parameter]WARNING! Deltat isn't set,  default value is taken as:"<< newmark.deltat << std::endl;
  }
  

  if (field_parser.search(3, "alpha","Alpha","ALPHA")) {
    newmark.alpha = field_parser.next(newmark.alpha);
  }
  else{
    newmark.alpha = 0;
    std::cout << "[Newmark Parameter]WARNING! Alpha isn't set,  default value is taken as:"<< newmark.alpha << std::endl;
  }
  

  if (field_parser.search(3, "beta","Beta","BETA")) {
    newmark.beta = field_parser.next(newmark.beta);
  }
  else{
    newmark.beta = 0.25;
    std::cout << "[Newmark Parameter]WARNING! Beta isn't set,  default value is taken as:"<< newmark.beta << std::endl;
  }
  
  if (field_parser.search(3, "gamma","Gamma","GAMMA")) {
    newmark.gamma = field_parser.next(newmark.gamma);
  }
  else{
    newmark.gamma = 0.5;
    std::cout << "[Newmark Parameter]WARNING! Gamma isn't set,  default value is taken as:"<< newmark.gamma << std::endl;
  }


  
};
}