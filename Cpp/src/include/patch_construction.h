/*
 * patch_construction.h
 *
 *  Created on: Apr 14, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_INTERSECTIONS_PARALLEL_PATCH_CONSTRUCTION_H_
#define COMMON_INTERSECTIONS_PARALLEL_PATCH_CONSTRUCTION_H_

#include "carl_headers.h"
#include "mesh_tables.h"

#include "intersection_tools.h"

namespace carl
{

// ************************
// Patch_Construction class
// ************************

/** \brief Class used to build a mesh patch from a parent mesh and an coupling mesh element.
 *  
 *    Other than saving and building the mesh patch, this class contains some access methods and
 *  data structures adapted to the advancing front intersection search method. These methods are marked with 
 *  \[ADV. FRONT\] (**Idea: transfer those to their own class ...**)
 *
 *    References:
 *
 *  \anchor Gander_article 1. Gander, M.J., Japhet, C.: An Algorithm for Non-
        Matching Grid Projections with Linear Complexity. In:
        M. Bercovier, M.J. Gander, R. Kornhuber, O. Widlund
        (eds.) Domain Decomposition Methods in Science and
        Engineering XVIII, no. 70 in Lecture Notes in Computational
        Science and Engineering, pp. 185{192. Springer
        Berlin Heidelberg (2009). URL http://link.springer.
        com/chapter/10.1007/978-3-642-02677-5_19. DOI:
        10.1007/978-3-642-02677-5 19
 */

class Patch_construction
{
protected:

  // Communicators and parallel data
  const libMesh::Parallel::Communicator&  m_comm;   ///< Parent mesh communicator.
  const unsigned int            m_nodes;  ///< Parent mesh number of processors.
  const unsigned int            m_rank;   ///< Parent mesh processor rank.

  const libMesh::Parallel::Communicator&  m_local_comm; ///< Patch mesh communicator

  // Meshes and point locators
  libMesh::ReplicatedMesh&          m_Mesh_parent;  ///< Parent mesh.
  libMesh::ReplicatedMesh           m_Mesh_patch; ///< Patch mesh.
  std::unique_ptr<libMesh::PointLocatorBase>  m_Parent_Point_Locator; ///< Parent mesh point locator.

  // Object used to do the intersection tests
  Intersection_Tools              m_Intersection_Test; ///< Intersection search and construction tools.

  // Data structures used to build the patch
  /// Set of elements inside the patch.
  std::unordered_set<unsigned int>                  m_Patch_Elem_indexes;
  /// Set of nodes inside the patch.
  std::unordered_set<unsigned int>                  m_Patch_Node_indexes;
  /// Element neighbour table for the patch.
  std::unordered_map<unsigned int,std::unordered_set<unsigned int> >  m_Patch_Elem_Neighbours;

  // Data structures used to do the advancing front intersection search (Adv. front)
  /// \[ADV. FRONT\] Index of the patch element currently being tested
  unsigned int              m_working_element_id; 

  /// !!!
  bool                  m_bTestNeighsForNewPairs;

  /** \brief \[ADV. FRONT\] Deque containing the elements to be treated.
   *
   *    In [\ref Gander_article "1"], this corresponds to either "bl" or "bil", depending on 
   *  whenever the patch is associated to the outer or inner 'while' loop.
   *
   */
  std::deque<int>             m_element_intersection_queue;

  /** \brief \[ADV. FRONT\] Deque containing the elements to be tested
   *
   *    In [\ref Gander_article "1"], this corresponds to "al", and is only used for the internal 'while' loop.
   *
   */
  std::deque<int>             m_element_test_queue;

  /** \brief \[ADV. FRONT\] Marks if an element was already treated or not
   *
   *    In [\ref Gander_article "1"], this corresponds to "bd" and "ad".
   *
   */
  std::unordered_map<unsigned int,int>    m_element_already_treated;

  /** \brief \[ADV. FRONT\] Marks if an element is already inside "m_element_intersection_queue"
   *
   *    This structure is used instead of a search inside "m_element_intersection_queue" due to 
   *  the inefficiency of doing a search inside a deque.
   */
  std::unordered_map<unsigned int,int>  m_element_inside_intersection_queue;

  /** \brief \[ADV. FRONT\] Set containing the current element's neighbors that must be tested yet.
   *
   *    In [\ref Gander_article "1"], this corresponds to "n".
   *
   */
  std::unordered_set<unsigned int>    m_element_neighbours_to_search;

  // Maps used to help building a patch mesh
  /// Mapping between the parent and patch mesh nodes, using the former as keys.
  std::unordered_map<unsigned int,int>        m_node_map_Parent_to_Patch;

  /// Mapping between the patch and parent mesh nodes, using the former as keys.
  std::unordered_map<unsigned int,int>        m_node_map_Patch_to_Parent;

  /// Mapping between the parent and patch mesh elements, using the former as keys.
  std::unordered_map<unsigned int,int>        m_elem_map_Parent_to_Patch;

  /// Mapping between the patch and parent mesh elements, using the former as keys.
  std::unordered_map<unsigned int,int>        m_elem_map_Patch_to_Parent;

  /// Print debug information? *Default:* false.
  bool m_bPrintDebug;

  // PROTECTED constructor
  Patch_construction(); //< Dummy constructor

  // PROTECTED methods
  /// Insert an parent mesh element inside the patch, updating the data structures.
  void insert_patch_element(const libMesh::Elem   * Patch_elem);

  /// Build a patch mesh from the patch data structures.
  void build_patch_mesh();

public:

  // Constructors
  /// Constructor with a pre-defined parent mesh and a local communicator.
  Patch_construction(libMesh::Mesh & mesh, const libMesh::Parallel::Communicator& local_comm, bool debugOutput = false) :
    m_comm { mesh.comm() },
    m_nodes { m_comm.size() },
    m_rank { m_comm.rank() },
    m_local_comm {  local_comm },
    m_Mesh_parent { mesh },
    m_Mesh_patch { m_local_comm },

    m_bPrintDebug { debugOutput }
  {
    m_Parent_Point_Locator = m_Mesh_parent.sub_point_locator();

    //  Instruction needed to avoid the code from crashing if a query is outside the mesh
    m_Parent_Point_Locator->enable_out_of_mesh_mode();

    m_working_element_id = 0;
    m_bTestNeighsForNewPairs = true;

    m_Patch_Elem_indexes.reserve(m_Mesh_parent.n_elem());
    m_Patch_Node_indexes.reserve(m_Mesh_parent.n_nodes());
  };

  // Getters
  libMesh::ReplicatedMesh & parent_mesh();      ///< Returns the parent mesh.
  libMesh::ReplicatedMesh & patch_mesh();       ///< Returns the patch mesh.
  std::unordered_set<unsigned int> & elem_indexes();  ///< Returns the set of elements inside the patch.
  std::unordered_set<unsigned int> & node_indexes();  ///< Returns the set of nodes inside the patch.
  unsigned int size();                ///< Returns the number of elements inside the patch.
  const libMesh::Elem * elem(unsigned int idx);   ///< Returns an element of the PARENT mesh.
  /// Returns the patch mesh element neighbour table.
  std::unordered_map<unsigned int,std::unordered_set<unsigned int> > & patch_elem_neighbours();

  // PUBLIC methods
  /** \brief Build the patch mesh covering a given "Query_elem".
   *  
   *    The patch is constructed using a variant of Gander's advancing front algorithm [\ref Gander_article "1"], using
   *  the "Query_elem" as the "mesh" of the outer 'while' loop, and the parent mesh in the inner 'while'
   *  loop. First, the list of elements from the parent mesh forming the patch is filled, and then the
   *  patch mesh is constructed using "Patch_construction::build_patch_mesh()"
   */
  void BuildPatch(const libMesh::Elem   * Query_elem);

  /// Convert an element index from the parent mesh to the patch mesh.
  unsigned int convert_parent_to_patch_elem_id(unsigned int input);

  /// Convert an element index from the patch mesh to the parent mesh.
  unsigned int convert_patch_to_parent_elem_id(unsigned int input);

  /// Export the patch mesh to a file.
  void export_patch_mesh(std::string & filename_base);

  unsigned int current_elem_id(); //< \[ADV. FRONT\] Get the current element index.
  const libMesh::Elem * current_elem_pointer(); //< \[ADV. FRONT\] Get the current element pointer.

  /*
   *    \[ADV. FRONT\] The following methods are all used by the advancing front algorithm,
   *  found inside the "intersection_search.h" files
   */

  /// \[ADV. FRONT\] Push back element to the deque containing the elements to be treated.
  void intersection_queue_push_back(unsigned int elem_id);

  /// \[ADV. FRONT\] Mark element as already treated.
  void set_elem_as_treated(unsigned int elem_id);

  /// \[ADV. FRONT\] Mark element as already inside the deque of elements to be treated.
  void set_elem_as_inside_queue(unsigned int elem_id);

  /// \[ADV. FRONT\] Check if deque of elements to be treated is empty.
  bool intersection_queue_empty();

  /// \[ADV. FRONT\] Check if deque of elements to be tested is empty.
  bool test_queue_empty();

  /// \[ADV. FRONT\] Pop and returns the first element to be treated.
  unsigned int intersection_queue_extract_front_elem();

  /// \[ADV. FRONT\] Pop and returns the first element to be tested.
  unsigned int test_queue_extract_front_elem();

  /// \[ADV. FRONT\] Returns the current element's neighbors that must be tested yet.
  std::unordered_set<unsigned int> & neighbors_to_search_next_pair();

  /** \brief \[ADV. FRONT\] Set the list of neighbors of the current element that must be tested yet.
   *
   *    Returns a boolean that short-circuits the search of new intersection pairs if the list is empty.
   */
  bool set_neighbors_to_search_next_pairs();

  /// \[ADV. FRONT\] Adds element to test list.
  void add_neighbors_to_test_list();

  /// \[ADV. FRONT\] Initialize the advancing front search data structures.
  void FrontSearch_initialize();

  /** \brief \[ADV. FRONT\] Reset the advancing front search data structures, with the exception of list of elements to be treated (**should be joined to "FrontSearch_prepare_for_probed_test"**)
   *
   *    This method is called before starting the inner loop of the advancing front algorithm.
   */
  void FrontSearch_reset();

  /** \brief \[ADV. FRONT\] Reset the advancing front search data structures, with the exception of list of elements to be treated.
   *
   *    This method is called before starting the inner loop of the advancing front algorithm.
   */
  unsigned int FrontSearch_prepare_for_probed_test();
};

}
#endif /* COMMON_INTERSECTIONS_PARALLEL_PATCH_CONSTRUCTION_H_ */
