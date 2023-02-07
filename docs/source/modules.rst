.. _modules:

Modules Description
###################

In this sections, we present the code structure, design and philosophy as well as the different modules. It does not details all the API and the different functions but the concepts and methods of the generic modules. The structure is based on 4 distinct types of module sometimes inter-connected: the ``connector``, the ``graph``, the ``wrap_helper`` and the ``algorithms``. 

These modules are connected via **standardised interfaces**, *i.e.* sets of predefined functions names and signatures that need to exists in any version of each module. 

The ``connector`` and the ``graph`` modules manage everything related to node connectivity, respectively at local and global level. All their functions helps user access to relationships between nodes and their geometry (e.g. getting all the downstream receivers,  getting the immediate neighbours, converting between node indices to X/Y coordinates, ...).


Name Convention and Indexing
============================

``DAGGER`` uses similar name conventions than `LANDLAB <https://landlab.readthedocs.io/en/master/user_guide/grid.html#basic-grid-elements>`_. **Nodes** are "points" connected to their neighbours via **links** and are at the centre of **cells** representing the node area of influence. Data, no matter what they represents, is either related to a node via a **node index** (e.g. elevation, total water flux) or a link via a **link index** (e.g. gradient, local water flux between two cells). In graph theory terms, nodes are the *vertices* and the links the *arcs*. One node typically has multiple links connecting to its neighbourhood.

Data is represented as 1D arrays, accessible *via* subsequent node indices or link indices.

TODO:Add a figure explaining that here

``Connector``
=============

The ``connector`` modules manages the topology of the grid at a node level. In other words, it manages the type of **grid** discretising the data and can be used in a stand-alone way. It deals with everything related to a given node and its immediate neighbouring: retrieving neighbour indices, related link indices, manages boundary conditions or no data, distance between two neighbours, local slope, partitioning weights ... Basically if you need any information concerning a node, its location and the relationship with its neighbours, this is the module to seek for. Anything "non local" will be managed by the ``graph``. Ultimately, the ``connector`` only needs geometrical plan-view information to fetch neighbours and links for each nodes (undirected graph): for example for regular grids it would be the number of row, col or the spacing in X and Y and the boundary codes (see bellow). If directionality is important, the ``connector`` can ingest a topographic field and compute a ``directed graph``, enabling the fetching of receivers and donors. Again the connector only bares this information at a node level - *i.e* which of the immediate neighbours are donors and receivers.

Boundary conditions
--------------------

``DAGGER`` defines a number of boundary condition types in ``DAGGER/boundary_conditions.hpp``. Before reading the details, note that default sets of boundary conditions are available in ``DAGGER`` and users only need to manually set them for specific cases. Boundary conditions are ``UINT8`` integers (under the form of an enumeration for the sake of clarity) for each and every nodes and are stored in an encapsulated class in the ``connector`` (``connector.boundaries``). The latter contains a 1D array of node size of boundary codes as well as helper functions. Essentially, they will define two separate aspects of boundaries: (i) whether a link can exist between two nodes and (ii) if flux can enter/leave a cell or just in one way. The boundary codes are:

- ``NO_FLOW (0)``: no data, the node will be ignored in most operations. Any link pointing to a ``NO_FLOW`` node is labelled "invalid", however ``NO_FLOW`` node can be listed when querying a list of neighbours. 
- ``FLOW (1)``: "normal" internal node. Flow can go through the cell in all directions but not leave the model (to the exception of explicitly unmanaged local minima).
- ``FLOW_BUT (2)``: Not used at the moment, but provided as an hypothetical boundary label for a node that can receive and give flux but user may want to differentiate it from a normal node to add additional checks. For example this could be used in the case of specific numerical treatment of node connected to a ``OUT`` boundary.
- ``CAN_OUT (3)``: nodes (often at the edge of the grid) where flux **can** leave the grid, but only if there is no downstream neighbours. 
- ``OUT (4)``: nodes (often at the edge of the grid) where flux **have to** leave the grid, even if there are downstream neighbours. Such nodes only receive flux but never give any.
- ``FORCE_OUT (5)``: nodes (often at the edge of the grid) where flux not only **have to** leave the grid, but they will force ANY neighbours that can give flux. Such nodes only receive flux but never give any.
- ``CANNOT_OUT (6)``: edge nodes where flux **CANNOT** leave the grid, they are basically ``FLOW`` nodes located at edges. The distinction from ``FLOW`` nodes remain important from a performance point of view: the code assumes a ``FLOW`` node has all its neighbours without checks while ``CANNOT_OUT`` nodes will spend extra computational effort to determine the exact number of neighbours (e.g. if the node is located at the first row of a regular grid).
- ``IN (7)``: nodes (often at the edge of the grid) that can only give flux to its downstream neighbours. They won't be able to receive flow in any cases.
- ``FORCE_IN (8)``: nodes (often at the edge of the grid) that not only give flux to its downstream neighbours, but also to ANY of their neighbour able to receive fluxes. They won't be able to receive flow in any cases.
- ``PERIODIC_BORDER (9)``: Edge nodes connected to the opposite boundary to simulate continuous flow lines. For example, in the case of a regular grid, a periodic node at the easternmost boundary would have neighbours on the Westernmost boundary.


``connector.boundaries`` class provide a number of helper function to check how a given node index has the possibility to manage flow (for example function ``is_normal_node`` is called like ``connector.is_normal_node(node_index)``):

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


Presets of boundary conditions are given (and detailed) for constructing the connector.

**IMPORTANT POINT**: these functions and boundary conditions defines the ability of a node to **potentially** create a link or receive/give/... data. Topography and boundary codes of neighbouring nodes determine if a link between the two is definitive.

