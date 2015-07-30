/*
 * common_header.h
 *
 *  Created on: Jul 29, 2015
 *      Author: Thiago Milanetto Schlittler
 *
 *      File containing all the external libraries includes (except the CGAL
 *	ones, included in "CGAL_typedefs")
 *
 */

#ifndef COMMON_HEADER_H_
#define COMMON_HEADER_H_

#include <unistd.h>

// --- Timing
#include <chrono>

// --- Boost/random
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/lagged_fibonacci.hpp>

// --- Boost/variant and functional
#include <boost/variant/apply_visitor.hpp>
#include <functional>

// --- IO
#include <fstream>
#include <iostream>

// --- Containers
#include <unordered_set>
#include <deque>
#include <vector>

// --- C++ strings
#include <string>

#endif /* COMMON_HEADER_H_ */
