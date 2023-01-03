Modules
===========

.. _modules:

In this sections, we present the code structure, design and philosophy as well as the different modules. It does not details all the API and the different functions but the concepts and methods of the generic modules. The structure is based on 4 distinct types of module sometimes inter-connected: the ``connector``, the ``graph``, the ``wrap_helper`` and the ``algorithms``. 

These modules are connected via **standardised interfaces**, *i.e.* sets of predefined functions names and signatures that need to exists in any version of each module. 

The ``connector`` and the ``graph`` modules manage everything related to node connectivity, respectively at local and global level. All their functions helps user access to relationships between nodes and their geometry (e.g. getting all the downstream receivers,  getting the immediate neighbours, converting between node indices to X/Y coordinates, ...).


Name Convention and Indexing
-----------------

``DAGGER`` uses similar name conventions than `LANDLAB <https://landlab.readthedocs.io/en/master/user_guide/grid.html#basic-grid-elements>`_. **Nodes** are connected to their neighbours via **links** and are at the centre of **cells**. Data, no matter what they represents, is either related to a node via a **node index** or a link via a **link index**. 

Data is represented as 1D arrays of number of node or number of link size. 

``Connector``
--------------

The ``Connector`` modules manages the topology of the grid at a node level. This module can be used in a stand-alone way and manages everything related to one node and its immediate neighbouring. For example this is where the number of neighbour for each node is defined, or the surface area each node represents, or the grid spacing to each link, ... 

It also manages the crucial aspect of boundary conditions: the ``connector`` can tell for each node whether "flow" can cross it, exit the model, be periodic or cyclic, ... 
