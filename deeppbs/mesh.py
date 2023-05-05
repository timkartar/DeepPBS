# Jared Sagendorf

"""A simple set of methods and classes for dealing with and performing basic operations on 
triangular meshes."""
# builtin modules
import os

# third party modules
import numpy as np
try:
    import igl
except BaseException as E:
    from .exceptions import ExceptionModule
    igl = ExceptionModule(E)

class Mesh(object):
    """Wrapper class for storing a trimesh mesh and peforming some basic operations"""
    def __init__(self, handle=None, vertices=None, faces=None, name="mesh", process=True, remove_disconnected_components=True, smoothing=None, smoothing_kwargs={}, **kwargs):
        try:
            import trimesh
        except ModuleNotFoundError:
            raise ModuleNotFoundError("Trimesh is required for mesh-related functionality!")
        
        self.name = name
        if handle is not None:
            if isinstance(handle, str):
                self.mesh = trimesh.load(handle, process=process, **kwargs)
            elif isinstance(handle, trimesh.Trimesh):
                self.mesh = handle
            else:
                raise TypeError("arg `handle` must be a string or a Trimesh object!")
        elif (vertices is not None) and (faces is not None):
            # Check that these are valid sizes/types
            if vertices.shape[1] != 3 or vertices.ndim != 2:
                raise ValueError("Vertices must be an Vx3 array.")
            if faces.shape[1] != 3 or faces.ndim != 2:
                raise ValueError("Faces must be an Fx3 array.")
            
            self.mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=process, **kwargs)
        else:
            raise Exception("Insuffucient data given to construct a mesh object.")
        
        if remove_disconnected_components:
            # remove any disconnected components
            self.remove_disconnected_components()
        
        if smoothing == "taubin":
            trimesh.smoothing.filter_taubin(self.mesh, **smoothing_kwargs)
        
        self.__reset()
    
    @property
    def area_faces(self):
        return self.mesh.area_faces
    
    @property
    def area(self):
        return self.mesh.area
    
    @property
    def volume(self):
        return self.mesh.volume

    @property
    def bbox(self):
        return self.mesh.bounding_box_oriented
    
    @property
    def aspect_ratio(self):
        return self.bbox.extents.max()/self.bbox.extents.min()
    
    @property
    def vertex_kdtree(self):
        return self.mesh.kdtree

    @property
    def vertex_adjacency_graph(self):
        return self.mesh.vertex_adjacency_graph
    
    @property
    def vertex_adjacency_matrix(self):
        if "vertex_adjacency_matrix" not in self.cache:
            try:
                from scipy.sparse import coo_matrix
                import networkx as nx
            except ModuleNotFoundError:
                raise ModuleNotFoundError("The dependencies 'networkx' and 'scipy' are required for this function!")
            
            self.cache["vertex_adjacency_matrix"] = coo_matrix(nx.adjacency_matrix(self.vertex_adjacency_graph, nodelist=[i for i in range(self.num_vertices)]))
        return self.cache["vertex_adjacency_matrix"]
    
    @property
    def undirected_edge_indices(self):
        """ Return a set of edge indices such that (i,j) and (j,i) are in E for
        all i and j pairs """
        if "undirected_edge_indices" not in self.cache:
            self.cache["undirected_edge_indices"] = np.stack([
                self.cache["vertex_adjacency_matrix"].row,
                self.cache["vertex_adjacency_matrix"].col
                ], axis=1
            )
        return self.cache["undirected_edge_indices"]
    
    @property
    def face_adjacency(self):
        return self.mesh.face_adjacency
    
    @property
    def vertex_normals(self):
        return self.mesh.vertex_normals
    
    @property
    def face_normals(self):
        return self.mesh.face_normals
    
    @property
    def vertex_faces(self):
        return self.mesh.vertex_faces
    
    @property
    def vertices(self):
        return self.mesh.vertices
    
    @property
    def faces(self):
        return self.mesh.faces
    
    @property
    def edges(self):
        return self.mesh.edges
    
    @property
    def convex_hull(self):
        return self.mesh.convex_hull
    
    @property
    def cot_matrix(self):
        if 'cot_matrix' not in self.cache:
            self.cache['cot_matrix'] = igl.cotmatrix(self.vertices, self.faces)
        return self.cache['cot_matrix']
    
    @property
    def mass_matrix(self):
        if 'mass_matrix' not in self.cache:
            self.cache['mass_matrix'] = igl.massmatrix(self.vertices, self.faces, igl.MASSMATRIX_TYPE_VORONOI)
        return self.cache['mass_matrix']
        
    @property
    def vertex_attributes(self):
        return self.mesh.vertex_attributes
    
    def __reset(self):
        """reset various mesh properties if something changes"""
        
        # simple properties of the mesh
        self.num_vertices = len(self.mesh.vertices) # number of vertices in the mesh
        self.num_faces = len(self.mesh.faces) # number of faces in the mesh
        self.num_edges = len(self.mesh.edges) # number of edges in the mesh
        
        # cached properties
        self.cache = {}
    
    def remove_disconnected_components(self):
        """Remove disconnected subcomponents keeping the largest"""

        components = self.mesh.split()
        sizes = [len(c.vertices) for c in components]
        
        self.mesh = components[np.argmax(sizes)]
        self.__reset()
    
    def nearestVertex(self, x):
        """Returns the distance and vertex index which is nearest the given point x"""
        d, i = self.vertex_kdtree.query(x)
        return i, d
    
    def verticesInBall(self, x, r):
        """Returns all verices within a radius 'r' of a point 'x'"""
        indices = self.vertex_kdtree.query_ball_point(x, r)
        indices = np.array(indices, dtype=np.int32)
        distances = np.linalg.norm(self.vertices[indices] - x, axis=1)
        
        return indices, distances
    
    def facesInBall(self, x, r):
        """Returns all triangles that contain at least one vertex within radius 'r' of point 'x'"""
        indices = self.vertex_kdtree.query_ball_point(x, r)
        fi = self.mesh.vertex_faces[indices].flatten()
        fi = np.unique(fi[fi >=0])
        
        return self.mesh.faces[fi]
    
    def findNeighbors(self, v, k=1, nlist=None, vo=None):
        """Returns the k-neighbors of a given vertex"""
        if nlist is None:
            nlist = set()
        if k == 0:
            return
        if vo is None:
            vo = v
        for n in self.vertex_adjacency_graph.neighbors(v):
            if n == vo:
                continue
            nlist.add(n)
            self.findNeighbors(n, k=k-1, nlist=nlist, vo=vo)
        return nlist
    
    def save(self, directory=".", file_name=None, file_format="off", overwrite=False):
        """Writes the mesh to file"""
        if file_name is None:
            file_name = "{}.{}".format(self.name, file_format)
        
        path = os.path.join(directory, file_name)
        exists = os.path.exists(path) # check if file already exists
        if exists and not overwrite:
            # just return the existing path and do not overwrite
            return path
        
        # save mesh to file
        self.mesh.export(file_name, file_format)

