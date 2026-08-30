#Cube file read and process
import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from itertools import combinations
from scipy.spatial import KDTree

class CubeFile:
    """
    Handles reading
    """
    BOHR2ANG = 0.529177

    def __init__(self, filename):
        self.filename = filename
        self.header = None
        self.natoms = None
        self.origin = None
        self.nval = 1
        self.n = []
        self.axis = []
        self.atoms = []
        self.elem = []
        self.data = None
        self._parse()
    
    def _parse(self):
        """Parse the cube file."""
        with open(self.filename,'r') as f:
            self.header = f.readline().split()
            line = f.readline().split()
            if line[0].isnumeric() == False:
                line = f.readline().split()
            self.natoms = int(line[0])
            self.origin = np.array(line[1:],dtype=float)*self.BOHR2ANG

            for i in range(3):
                line = f.readline().split()
                self.n.append(int(line[0]))
                self.axis.append(np.array(line[1:],dtype=float)*self.BOHR2ANG)

            for i in range(self.natoms):
                line = f.readline().split()
                self.elem.append(int(line[0]))
                self.atoms.append(np.array(line[2:],dtype=float)*self.BOHR2ANG)

            remaining = f.read()
            self.data = np.fromstring(remaining, dtype=float, sep=' ').reshape(
                self.n[1],self.n[0],self.n[2]).transpose(1,0,2)   

    
    def get_coordinates(self,fractional=False):
        """Generate physical coordinates for all grid points."""
        grid_coords = []
        shape = self.n
        axis = self.axis
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    # Coordinate of the voxel corner
                    if fractional==False:
                        coord = self.origin + (i * axis[0] +
                                            j * axis[1] + k * axis[2])
                    else:
                        coord = [i/(shape[0]),j/(shape[1]),k/(shape[2])]
                    weight = self.data[i, j, k] # type: ignore
                    grid_coords.append([coord[0], coord[1], coord[2], weight])
        return np.array(grid_coords)

    def frac_to_absolute(self,coord):
            
        return
    
class ProcessGrid():
    '''
    Process the Cube grid in xyzw format
    '''
    def __init__(self, cube:CubeFile):
        self.data = cube.get_coordinates()
        self.cube = cube.data
        self.axis = np.array(cube.axis)
        self.n = cube.n
        self.axis_scaled = self.axis * self.n
        self.frac_coord = cube.get_coordinates(fractional=True)
        self.max_frac = np.max(self.frac_coord[:,:-1],axis=0)
        self.interpolator = None
        
        # Pre-compute absolute grid positions
        self.grid = self.data[:, :-1]

    def Nearest_Coord(self,coord):
        subtract = []
        for i in coord:
            if i > 0:
                subtract.append(np.ceil(i))
            else:
                subtract.append(np.floor(i))
        subtract = np.array(subtract)

        return self.frac_coord[np.argmin(np.linalg.norm(self.frac_coord[:,:-1]+subtract-coord,axis=1)),:-1] + subtract

    @staticmethod
    def frac_to_abs(coord:np.ndarray,axes):
        if len(coord.shape) == 1 or coord.shape[0] == 1:
            return coord[0]*axes[0] + coord[1]*axes[1] + coord[2] * axes[2]
        else:
            return np.outer(coord[:, 0], axes[0]) + np.outer(coord[:, 1], axes[1]) + np.outer(coord[:, 2], axes[2])

    @staticmethod
    def abs_to_frac(coord,axes):
        coord = np.array(coord)
        return (np.linalg.inv(axes) @ coord.T).T

    def map_to_cell(self,coord):
        wrap = lambda x, tol=1e-12: np.where((np.abs(x % 1.0) < tol) | (np.abs(x % 1.0 - 1.0) < tol), 0.0,
        x % 1.0)
        coord_wrapped = wrap(self.abs_to_frac(coord,self.axis_scaled))
        coord_wrapped = self.frac_to_abs(coord_wrapped,self.axis_scaled)
        dif = coord_wrapped - coord
        _, indices = self.kdtree.query(coord_wrapped)
        closest = self.data[indices] 
        closest[:,:-1] = closest[:,:-1] - dif
        return closest

    @staticmethod
    def one_arr(i):
        i = int(i)
        arr = np.array([0,0,0])
        arr[i] = 1
        return arr
    
    def Expand_abs(self,min,max=[]):
        if len(max) == 0:
            max = self.max_frac
        else:
            max = np.array(max)
        origin_frac = self.Nearest_Coord(np.array(min))
        dist = np.array((max - origin_frac)*self.n[0],dtype=int)

        origin_abs = self.frac_to_abs(origin_frac,self.axis)
        new_grid = []

        wrap = lambda x, tol=1e-12: 0.0 if abs(x % 1.0) < tol or abs(x % 1.0 - 1.0) < tol else x % 1.0
        for i in range(dist[0]):
                    for j in range(dist[1]):
                        for k in range(dist[2]):
                            # Coordinate of the voxel corner+
                            coord = origin_frac +np.array([i/(self.n[0]),j/(self.n[1]),k/(self.n[2])])
                            weight = self.cube[int(wrap(coord[0])*self.n[0]), int(wrap(coord[1])*self.n[1]),int(wrap(coord[2])*self.n[2]) ] # type: ignore
                            coord = self.frac_to_abs(coord,self.axis)
                            new_grid.append([coord[0], coord[1], coord[2], weight]) # type: ignore
        return new_grid
    
    def Expand_frac(self,min,max=[]):
        if max == []:
            max = self.max_frac
        else:
            max = np.array(max)
        origin_frac = self.Nearest_Coord(np.array(min))
        dist = np.array((max - origin_frac)*self.n[0],dtype=int)
        new_grid = []
        for i in range(dist[0]):
                    for j in range(dist[1]):
                        for k in range(dist[2]):
                            # Coordinate of the voxel corner
                            
                            coord = origin_frac + np.array([i/(self.n[0]),j/(self.n[1]),k/(self.n[2])])
                            weight = self.cube[int(coord[0]%1)*self.n[0], int(coord[1] %1)*self.n[1], int(coord[2]%1)*self.n[2] ] # type: ignore
                            new_grid.append([coord[0], coord[1], coord[2], weight])
        return new_grid
    

    def around_point(self,coord,length=1):
        # Pre-compute inverse matrix once
        self.inv_axes = np.linalg.inv(self.axis_scaled)
        
        # Pre-compute wrapped grid points in fractional coordinates
        # This avoids repeated wrapping of the grid
        self.grid_frac = self.abs_to_frac(self.grid,self.axis_scaled)
        
        # Build KD-tree for fast nearest neighbor search
        self.kdtree = KDTree(self.grid)

        npt = 50 
        while True:
            x = np.linspace(coord[0]-length/2,coord[0]+length/2, npt)
            y = np.linspace(coord[1]-length/2,coord[1]+length/2, npt)
            z = np.linspace(coord[2]-length/2,coord[2]+length/2, npt)
            # Create meshgrid
            X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

            # Combine into array of [x, y, z] points
            combined = np.stack([X, Y, Z], axis=-1).reshape(-1, 3)
            mapped = self.map_to_cell(combined)
            if mapped.shape != np.unique(mapped,axis=0).shape:
                npt -= 1
                print(npt)
            else:
                print(npt,'OK')
                return mapped
        
c = CubeFile("defc_afm1_up.cube")
print(c.atoms[2])
grid = ProcessGrid(c)
g = grid.around_point(c.atoms[2],length=3) # type: ignore
np.savetxt('test1.txt',g)
# g=grid.Expand_abs(min=[-0.1,0,0])
pass

