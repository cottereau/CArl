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

// --- Timing
#include <chrono>

// --- Boost/random
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/lagged_fibonacci.hpp>

// --- Boost/variant and functional
//#include <boost/variant/apply_visitor.hpp>
//#include <functional>

// --- Boost/hash
//#include <boost/functional/hash.hpp>

// --- Boost/filesystem
//#include <boost/filesystem.hpp>

// --- IO
#include <fstream>
#include <iostream>

// --- Containers and algorithms
#include <unordered_set>
#include <unordered_map>
#include <deque>
#include <vector>
#include <iterator>
#include <algorithm>
#include <tuple>
#include <map>

// --- C++ strings
#include <string>
#include <sstream>

// --- Common C/C++
#include <math.h>
#include <unistd.h>

#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <array>

extern boost::random::lagged_fibonacci607 m_rng;

// A small homemade assert
#define homemade_assert_msg(asserted, msg) do \
{ \
	if (!(asserted)) \
	{ \
		std::cerr << "Assertion `" #asserted "` failed in " << __FILE__ << std::endl\
				<< " > line " << __LINE__ << ": " <<  msg << std::endl; \
		std::exit(EXIT_FAILURE); \
	} \
} while(false)

#define homemade_error_msg(msg) do \
{ \
	std::cerr << "Error: " <<  msg << std::endl; \
	std::exit(EXIT_FAILURE); \
} while(false)

#endif /* COMMON_HEADER_H_ */
