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

// *************************
// Intersection_Search class
// *************************

/** \brief Class containing the structure needed to find all the
 *  intersections between two meshes, inside the coupling region mesh.
 */

class Intersection_Search
{
protected:

  // Mesh addresses
  libMesh::Mesh&           m_Mesh_A;
  libMesh::Mesh&           m_Mesh_B;
  libMesh::Mesh&           m_Mesh_Coupling; ///< Mesh representing the coupling region

  // Communicators and parallel data
  const libMesh::Parallel::Communicator& m_comm;      ///< Global communicator
  const unsigned int             m_nodes;     ///< Number of MPI processors
  const unsigned int             m_rank;      ///< This processor's rank (or ID in the communicator)
  const libMesh::Parallel::Communicator& m_local_comm;  ///< Local communicator, used for mesh patches

  // Objects containing the two patches
  Patch_construction             m_Patch_Constructor_A; ///< Patch_construction object for mesh A
  Patch_construction             m_Patch_Constructor_B; ///< Patch_construction object for mesh B

  /** \brief Flag defining the intersection meshing algorithm.
   *
   * Can be either carl::LIBMESH_TETGEN (use LibMesh's tetgen module, problematic with Intel compilers) 
   * or carl::CGAL (use CGAL's Triangulation_3). *Default:* carl::CGAL.
   */
  IntersectionMeshingMethod m_MeshingMethod;  

  /// \brief Object containing the intersection mesh.
  Mesh_Intersection            m_Mesh_Intersection;

  /** \brief Multimap containing the intersection pairs.
   *
   * The element indexes of the mesh A are used as the keys of the multimap
   */
  std::unordered_multimap<unsigned int,unsigned int> m_Intersection_Pairs_multimap;

  /// \brief Intersection operations for the main intersection tests
  Intersection_Tools m_Intersection_test;

  /// \brief Intersection operations for the neighbor intersection tests. (**why do we need both?**)
  Intersection_Tools m_Intersection_test_neighbors;

  /// \brief Vector containing the number of intersections found inside each of the coupling mesh elements.
  std::vector<unsigned int> m_Nb_Of_Intersections_Elem_C;

  /// \brief libMesh::ErrorVector used to repartition the coupling mesh.
  libMesh::ErrorVector m_coupling_weights;

  /// \brief Boolean flag determining if we should save the intersection data or not. *Default:* true.
  bool m_bSaveInterData;

  /// \brief Boolean flag determining if we did a preliminary intersection search, used for the coupling mesh repartition and memory preallocations.
  bool m_bDidPreliminarySearch;

  /// \brief Boolean flag determining if we have the data for a proper preallocation of m_Intersection_Pairs_multimap.
  bool m_bHavePreallocData;

  /// \brief Boolean flag indicating if the intersections were built. (**A bit useless right now ...**)
  bool m_bIntersectionsBuilt;

  /// \brief Volume cutoff for the intersection polyhedrons. *Default:* 1E-15.
  double m_Min_Inter_Volume;

  /// \brief Boolean indicating if we should skip the intersection construction. *Default:* false.
  bool m_bSkipIntersectionConstruction;

  /// \brief Boolean indicating if we should skip the intersection repartitioning. *Default:* false.
  bool m_bSkipIntersectionPartitioning;

  std::string m_Output_filename_base;     ///< Output filenames base.

  bool MASTER_bPerfLog_intersection_search; ///< Do performance log? *Default:* true.
  libMesh::PerfLog m_perf_log;        ///< libMesh::PerfLog object.
  bool m_bPrintDebug;             ///< Print debug data. *Default:* false.
  bool m_bPrintTimingData;          ///< Print timing data. *Default:* false.
  bool m_bPrintIntersectionsPerPartData;    ///< Print intersections per partition. *Default:* false.
  std::string m_timing_data_file_base;    ///< Output filenames base for the timing data. 

  // PROTECTED methods
  /** \brief Build both patches associated a given query element from the coupling region mesh.
   *    
   * Calls the method Patch_construction::BuildPatch for both patches, and export the patch meshes if Intersection_Search::m_bPrintDebug == true.
   */ 
  void BuildCoupledPatches(const libMesh::Elem  * Query_elem, int patch_counter = 0);

  /** \brief Find all the intersections between the patches associated to the Query_elem, 
   * using a brute force method. All elements from a patch are tested against all the elements from the other 
   * patch.
   */
  void FindPatchIntersections_Brute(const libMesh::Elem   * Query_elem);

  /** \brief Find all the intersections between the patches associated to the Query_elem, using an advancing
   *  front method.
   */
  void FindPatchIntersections_Front(const libMesh::Elem   * Query_elem);

  /** \brief Find the first intersecting pair from a patch.
   *  
     *  Do so by :
   *    1. locating the elements from the guide patch that contain the
   *    vertices of an element from the probed patch.
   *    2. for each one of these guide elements, find the one that
   *    intersects the probed element inside the coupling.
   *    3. if this fails, use a brute force search method (call Intersection_Search::BruteForce_FindFirstPair).
   */
  void FindFirstPair(   Patch_construction *  Patch_guide,
              Patch_construction *  Patch_probed,
              std::pair<unsigned int,unsigned int> &    First_intersection);

  /** \brief Find the first intersecting pair from a patch, doing a full scan of
   *  the patches.
   *
   *  This is, essentially, a brute force algorithm set to stop after the
   *  first positive test).
   */
  void BruteForce_FindFirstPair(  Patch_construction *  Patch_guide,
                  Patch_construction *  Patch_probed,
                  std::pair<unsigned int,unsigned int> &    First_intersection);

  /** \brief Find and build all the intersections, using the brute force method.
   *
   * Skips build step if Intersection_Search::m_bSkipIntersectionConstruction == true.
   */
  void FindAndBuildIntersections_Brute();

  /** \brief Find all the intersections, using the brute force method.
   *
   * Skips saving intersecting element pairs if Intersection_Search::m_bSaveInterData == true.
   */
  void FindIntersections_Brute();

  /** \brief Find and build all the intersections, using the advancing front method.
   *
   * Skips build step if Intersection_Search::m_bSkipIntersectionConstruction == true.
   */
  void FindAndBuildIntersections_Front();

  /** \brief Find all the intersections, using the advancing front method.
   *
   * Skips saving intersecting element pairs if Intersection_Search::m_bSaveInterData == true.
   */
  void FindIntersections_Front();

  /// \brief Calculate the volume of the intersections associated to Query_elem without updating the intersection mesh. **Currently unused**
  void CalculateIntersectionVolume(const libMesh::Elem  * Query_elem);

  /// \brief Build the intersections associated to the Query_elem and update the intersection mesh.
  void UpdateCouplingIntersection(const libMesh::Elem   * Query_elem);

public:

  /// \brief Constructor.
  Intersection_Search(libMesh::Mesh & mesh_A,
            libMesh::Mesh & mesh_B,
            libMesh::Mesh & mesh_Coupling,
            libMesh::Mesh & mesh_I,
            const std::string & output_base = std::string("test"),
            IntersectionMeshingMethod MeshingMethod = IntersectionMeshingMethod::CGAL,
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
    m_Patch_Constructor_A { m_Mesh_A,m_local_comm},
    m_Patch_Constructor_B { m_Mesh_B,m_local_comm},
    m_MeshingMethod { MeshingMethod },
    m_Mesh_Intersection { Mesh_Intersection(mesh_I,m_Mesh_A,m_Mesh_B,m_MeshingMethod)},
    m_bSaveInterData { true },
    m_bDidPreliminarySearch { false },
    m_bHavePreallocData { false },
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
  };

  /// \brief Reference to the mesh A.
  libMesh::Mesh & mesh_A();

  /// \brief Reference to the mesh B.
  libMesh::Mesh & mesh_B();

  /// \brief Reference to the coupling mesh.
  libMesh::Mesh & mesh_Coupling();

  // PUBLIC Methods
  /** \brief Do a preliminary search, to optimize the intersection search and construction.
   *
   * This method finds all the intersecting element pairs, only saving their number for each
   * element of the coupling mesh. This information is used to repartition the coupling mesh 
   * and to preallocate the Intersection_Search::m_Intersection_Pairs_multimap.
   */
  void PreparePreallocationAndLoad(SearchMethod search_type = BRUTE);

  /// \brief Preallocate Intersection_Search::m_Intersection_Pairs_multimap and repartition the coupling mesh.
  void PreallocateAndPartitionCoupling();

  /// \brief Find and build all the intersections.
  void BuildIntersections(SearchMethod search_type = BRUTE);

  /// \brief Set up timing output file.
  void SetScalingFiles(const std::string& timing_data_file_base);

  /// \brief Set the "Skip intersection construction" flag.
  void SkipIntersectionConstruction(bool bSkipIntersectionConstruction)
  {
    m_bSkipIntersectionConstruction = bSkipIntersectionConstruction;
  }

  /// \brief Set the "Skip the coupling mesh repartition" flag.
  void SkipIntersectionPartitioning(bool bSkipIntersectionPartitioning)
  {
    m_bSkipIntersectionPartitioning = bSkipIntersectionPartitioning;
  }

  /// \brief Calculate the total intersection volume over all the processors
  void CalculateGlobalVolume();
};
}

#endif /* COMMON_INTERSECTIONS_PARALLEL_INTERSECTION_SEARCH_H_ */
