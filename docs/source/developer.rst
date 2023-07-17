.. _developer:

Developer guide
===============


|You can find here information for extending or helping to the development of ``DAGGER``. The global philosophy is straightforward: each module follow a **standard interface**, i.e. a list of attribute and methods the given module **needs** to be compatible with the core of the framework.

|Of course, they can have as many additional functionalities as wanted - the standard interfaces are kept minimal to enforce unconditional compatibility between modules. Think of it like a car dashboard and its engine. No matter how fancy or cheap your dashboard is, there needs to be the same connectivy with the engine. THe rest of the dashboard are just optional.

.. Note:: ``Connector``s need to be ``c++`` classes, where  **attributes** refer to class member variables and **methods** refer to class member function


Standard Interfaces for the Connectors
---------------------------------------

In order to be compatible with the ``graph``, a ``connector`` needs the following **attributes** (type in brackets, vector are ``std:vector`` of [type, size]):

- ``nnodes`` (int): the total number of nodes
- ``nnodes_t`` (size_t): the total number of nodes as unsigned integer
- ``nneighbours`` (int): the maximum number of neighbour for one node
- ``nneighbours_t`` (size_t): the maximum number of neighbour for one node as unsigned integer
- ``Sreceivers`` (vector [int, nnodes]): index is node index and value is the node index of the steepest descent receiver (once computed).
- ``Sdistance2receivers`` (vector [float, nnodes]): index is node index and value is the distance to the steepest descent receiver (once computed).
- ``Sndonors`` (vector [int, nnodes]): index is node index and value is the number of steepest descent donors, *i.e. the number of nodes having the given node index as Sreceivers* (once computed).
- ``links`` (vector [uint8, nlinks]): vector of link type (see nodule documentation for meaning).
- ``linknodes`` (vector [int, nlinks * 2]): vector of link nodes where link index * 2 and linkindex * 2 + 1 are the nodes composing a link.
- ``boundaries`` (boundary object): see boundary standards


It also need the following **methods** (written in ``c++`` format: return type neameoffunction(param1, param2)):

- ``void _allocate_vectors()``: initialise data structure each time the size of the grid changes (allocate memory for all the vectors needing it)
- ``void _reallocate_vectors()``: reinitialise vectors to default value
- ``void fill_linknodes()``: Uses the boundary conditions to (re)initialise the ``linknodes`` arrays, if explicitly stored by the ``connector`` and not recomputed at each requests.
- ``void update_links(array1D topography)``: uses a topographic array to determine link types (normal, inverse, invalid) accordingly and calculate steepest descent receivers/donors/...
- ``void update_links_MFD_only(array1D topography)``: uses a topographic array to determine link types (normal, inverse, invalid) accordingly. Does not compute any steepest descent info.
- ``void recompute_SF_donors_from_receivers()``: if steepest receivers donors are explicitely stored in the graph, recompute the list based on the ``Sreceivers`` array
- ``void get_nth_steepest_donors_of_node(int node, int nneighbours)``: returns the nth steepest descent donors of a given node. Might sound oddly specific, but is a crucial component of the very efficient SFD topological sorting (modified from Braun and Willett, 2013)
- ``vector[int, nneighbours] get_empty_neighbour()``: returns an empty vector of max nneighbours size, used to initialised neighbours/receivers/donors fetchers.
- ``int get_receivers_idx(int node, vector[int, nneihbours]& recs)``: take a node index, fetch its receivers (precomputed by an update_links function!), store them **in place** in the recs array and returns the number of receivers. While not the most instinctive function, it is by far the most efficient way to repeatedly fetch receivers info.
- ``int get_donors_idx(int node, vector[int, nneihbours]& dons)``: Same idea, but with donors.
- ``int get_neighbours_idx(int node, vector[int, nneihbours]& dons)``: Same idea, but with all neighbours.
- ``bool flow_out_model(int node)``: uses the links and boundary condition information to determine whether flow actually leaves of the model or not.
- ``bool flow_out_or_pit(int node)``: Same idea, but also returns true if the node is a pit.
- ``bool is_link_valid(int link_index)``: Check if a given link is valid or not. As link index is usually linked to node index via a geometrical relationship, some link indices can theoretically exist while not being valid depending on boundary conditions.
- ``int get_from_link(int link_index)``: take a link index and returns the node on the donor side.
- ``int get_to_link(int link_index)``: take a link index and returns the node on the receiver side.
- ``bool is_in_bound(int node_index)``: take a node index and check whether it exists or not.
- ``float get_dx_from_links_idx(int link_index)``: returns the distance represented by a link from its link id


Standard Interfaces for the Boundary object
--------------------------------------------

The boundary condition object is tightly connected to the ``connector``. Its format is flexible, as long as it can answer the following queries:

- ``is_normal_node``: returns true if the boundary is ``FLOW``
- ``is_bc``:returns true if the node needs specific boundary care for computing neighbours (i.e. is not ``FLOW`` or ``FLOW_BUT``)
- ``can_receive``: returns true if flux can enter a node from another one
- ``can_give``: returns true if flux can leave a node to another one
- ``can_flow_through``: returns true if a node can receive and give flow
- ``can_out``: returns true if the flux can leave the model *via* this cell
- ``no_data``: returns true if the node is ``NO_FLOW``
- ``can_create_link``: returns true if the node can **POTENTIALLY** create links with neighbours (it does not mean the link will be created, it may depend on the topography and/or the boundary code of the other node - which is not known by this function)
- ``forcing_io``: returns true if this nodes forces a link direction no matter the topography (e.g. ``FORCE_IN, FORCE_OUT``)
- ``force_giving``: returns true if this nodes forces flow to neighbours (e.g. ``FORCE_IN``)
- ``force_receiving``: returns true if this nodes forces flow from neighbours (e.g. ``FORCE_OUT``)
- ``is_periodic``: is the boundary code forcing periodicity
- ``normal_neighbouring_at_boundary``: returns true if the neighbouring is "normal" for a node (as opposed to periodic)
