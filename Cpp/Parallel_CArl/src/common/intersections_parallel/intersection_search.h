/*
 * intersection_search.h
 *
 *  Created on: Apr 14, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_SEARCH_H_
#define COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_SEARCH_H_

#include "carl_headers.h"
#include "mesh_tables.h"

#include "mesh_intersection_methods.h"
#include "patch_construction.h"

namespace carl
{

enum SearchMethod
{
	BRUTE,
	FRONT,
	BOTH
};

/*
 * 		Intersection_Search class
 *
 * 			This class contains the structure needed to find all the
 * 		intersections between two meshes A and B, and the coupling region mesh
 * 		(C or Coupling).
 *
 */

class	Intersection_Search
{
protected:

	// Meshes
	libMesh::Mesh&				   m_Mesh_A;
	libMesh::Mesh&				   m_Mesh_B;
	libMesh::Mesh&				   m_Mesh_Coupling;

	// Communicators and parallel data
	/*
	 * 		The local variants are used for the patch meshes
	 */
	const libMesh::Parallel::Communicator& m_comm;
	const unsigned int					   m_nodes;
	const unsigned int					   m_rank;
	const libMesh::Parallel::Communicator& m_local_comm;

	// Objects containing the two patches
	Patch_construction					   m_Patch_Constructor_A;
	Patch_construction					   m_Patch_Constructor_B;

	// Object containing the intersection mesh
	Mesh_Intersection					   m_Mesh_Intersection;

	// Output multimap containing the intersection pairs
	std::unordered_multimap<unsigned int,unsigned int> m_Intersection_Pairs_multimap;

	// Objects used to test the intersections
	Intersection_Tools m_Intersection_test;
	Intersection_Tools m_Intersection_test_neighbors;

	// Vector saving the number of intersections found inside each of the
	// coupling mesh elements
	std::vector<unsigned int> m_Nb_Of_Intersections_Elem_C;
	libMesh::ErrorVector m_coupling_weights;

	// Boolean flag determining if we should save the intersection data or not
	bool m_bSaveInterData;

	// Boolean flag determining if we did a preallocation run or not
	bool m_bPreparedPreallocation;

	// Boolean flag determining if we did a preallocation run or not
	bool m_bDidPreallocation;

	// Boolean flag indicating if the intersections were built
	bool m_bIntersectionsBuilt;

	// Volume cutoff for the intersections
	double m_Min_Inter_Volume;

	// Boolean indicating if intersections must be built or not
	bool m_bSkipIntersectionConstruction;
	bool m_bSkipIntersectionPartitioning;

	// Output
	std::string m_Output_filename_base;

	// Perflog and debug variables
	bool MASTER_bPerfLog_intersection_search;
	libMesh::PerfLog m_perf_log;
	bool m_bPrintDebug;
	bool m_bPrintTimingData;
	bool m_bPrintIntersectionsPerPartData;
	std::string m_timing_data_file_base;

	// PROTECTED methods
	/*
	 * 		Build both patches associated to the query element
	 */
	void BuildCoupledPatches(const libMesh::Elem 	* Query_elem, int patch_counter = 0);

	/*
	 * 		Find all the intersections between the patches, using a brute force
	 * 	method (all elements from a patch are tested against all the elements
	 * 	from the other patch).
	 */
	void FindPatchIntersections_Brute(const libMesh::Elem 	* Query_elem);

	/*
	 * 		Find the first intersecting pair from a patch. Do so by :
	 *
	 * 		a) locating the elements from the guide patch that contain the
	 * 		vertices of an element from the probed patch.
	 * 		b) for each one of these guide elements, find the one that
	 * 		intersects the probed element inside the coupling.
	 * 		c) if this fails, use a brute force search method.
	 */
	void FindFirstPair(		Patch_construction * 	Patch_guide,
							Patch_construction * 	Patch_probed,
							std::pair<unsigned int,unsigned int> &		First_intersection);

	/*
	 * 		Find the first intersecting pair from a patch, doing a full scan of
	 * 	the patches (essentially, a brute force algorithm set to stop after the
	 * 	first positive test).
	 */
	void BruteForce_FindFirstPair(	Patch_construction * 	Patch_guide,
									Patch_construction * 	Patch_probed,
									std::pair<unsigned int,unsigned int> &		First_intersection);

	/*
	 * 		Find all the intersections between the patches, using an advancing
	 * 	front method.
	 */
	void FindPatchIntersections_Front(const libMesh::Elem 	* Query_elem);

	/*
	 * 		For each coupling element, build the patches and find their
	 * 	intersections, using the brute force method.
	 */
	void BuildIntersections_Brute();
	void PrepareIntersections_Brute();

	/*
	 * 		For each coupling element, build the patches and find their
	 * 	intersections, using the advancing front method.
	 */
	void BuildIntersections_Front();
	void PrepareIntersections_Front();

	/*
	 * 		Legacy function, used to calculate the volume of the intersections
	 * 	without updating the intersection mesh.
	 */
	void CalculateIntersectionVolume(const libMesh::Elem 	* Query_elem);

	/*
	 * 		Take the intersection tables info and update the intersection mesh.
	 */
	void UpdateCouplingIntersection(const libMesh::Elem 	* Query_elem);

public:

	// Constructor
	Intersection_Search(libMesh::Mesh & mesh_A,
						libMesh::Mesh & mesh_B,
						libMesh::Mesh & mesh_Coupling,
						libMesh::Mesh & mesh_I,
						const std::string & output_base = std::string("test"),
						double Min_Inter_Volume = 1E-15,
						bool  bDoPerf_log = true,
						bool  bDebugOutput = false) :
		m_Mesh_A { mesh_A },
		m_Mesh_B { mesh_B },
		m_Mesh_Coupling { mesh_Coupling },
		m_comm { m_Mesh_Coupling.comm() },
		m_nodes { m_comm.size() },
		m_rank { m_comm.rank() },
		m_local_comm { mesh_I.comm() },
		m_Patch_Constructor_A { Patch_construction(m_Mesh_A,m_local_comm)},
		m_Patch_Constructor_B { Patch_construction(m_Mesh_B,m_local_comm)},
		m_Mesh_Intersection { Mesh_Intersection(mesh_I,m_Mesh_A,m_Mesh_B)},
		m_bSaveInterData { true },
		m_bPreparedPreallocation { false },
		m_bDidPreallocation { false },
		m_bIntersectionsBuilt { false },
		m_Min_Inter_Volume { Min_Inter_Volume },
		m_bSkipIntersectionConstruction { false },
		m_bSkipIntersectionPartitioning { false },
		m_Output_filename_base { output_base + "_r_" + std::to_string(m_rank) + "_n_" + std::to_string(m_nodes)},
		MASTER_bPerfLog_intersection_search {bDoPerf_log},
		m_perf_log { libMesh::PerfLog("Intersection search", MASTER_bPerfLog_intersection_search) },
		m_bPrintDebug { bDebugOutput },
		m_bPrintTimingData { false },
		m_bPrintIntersectionsPerPartData { false }
	{
		// Reserve space for the intersection multimap
		m_Nb_Of_Intersections_Elem_C.resize(mesh_Coupling.n_elem(),0);
		// m_Intersection_Pairs_multimap.reserve(mesh_A.n_elem()*mesh_B.n_elem());
	};

	// Getters
	libMesh::Mesh & mesh_A();
	libMesh::Mesh & mesh_B();
	libMesh::Mesh & mesh_Coupling();

	// PUBLIC Methods
	/*
	 * 		Preallocate run. It essentially does the intersection run, but
	 * 	without saving the data or building the intersections themselves.
	 *
	 */
	void PreparePreallocationAndLoad(SearchMethod search_type = BRUTE);
	void PreallocateAndPartitionCoupling();

	/*
	 * 		Interface for the user to build the intersections. By default, it
	 * 	uses the brute force algorithm, but the argument can be changed to
	 * 	carl::FRONT to use the advancing front method, of to carl::BOTH to use
	 * 	both methods (useful for benchmarking).
	 */
	void BuildIntersections(SearchMethod search_type = BRUTE);

	/*
	 * 		Set up timing output file.
	 */
	void SetScalingFiles(const std::string& timing_data_file_base);

	/*
	 * 		Skip intersection construction (or not)
	 */
	void SkipIntersectionConstruction(bool bSkipIntersectionConstruction)
	{
		m_bSkipIntersectionConstruction = bSkipIntersectionConstruction;
	}

	void SkipIntersectionPartitioning(bool bSkipIntersectionPartitioning)
	{
		m_bSkipIntersectionPartitioning = bSkipIntersectionPartitioning;
	}

	/*
	 * 		Calculate the volume over all the processors
	 */
	void CalculateGlobalVolume();
};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_SEARCH_H_ */
