#Cube file read and process
import numpy as np

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
                self.atoms.append(np.array(line[1:],dtype=float))

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
                        coord = self.origin + (i * axis[0]/shape[0] +
                                            j * axis[1]/shape[1] + k * axis[2]/shape[2])
                    else:
                        coord = [i/(shape[0]),j/(shape[1]),k/(shape[2])]
                    weight = self.data[i, j, k]
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
        self.axis = cube.axis
        self.n = cube.n
        self.frac_coord = cube.get_coordinates(fractional=True)
        self.max_frac = np.max(self.frac_coord[:,:-1],axis=0)

    def Nearest_Coord(self,coord):
        subtract = []
        for i in coord:
            if i > 0:
                subtract.append(np.ceil(i))
            else:
                subtract.append(np.floor(i))
        subtract = np.array(subtract)

        return self.frac_coord[np.argmin(np.linalg.norm(self.frac_coord[:,:-1]+subtract-coord,axis=1)),:-1] + subtract

    def frac_to_abs(coord:np.ndarray,axes:np.ndarray):
        if len(coord.shape) == 1 or coord.shape[0] == 1:
            return coord[0]*axes[0] + coord[1]*axes[1] + coord[2] * axes[2]
        else:
            return coord[:,0]*axes[0] + coord[:,1]*axes[1] + coord[:,2]*axes[2]
        
    def Expand_abs(self,min,max=[]):
        if max == []:
            max = self.max_frac
        else:
            max = np.array(max)
        max = 
        origin_new = self.Nearest_Coord(np.array(min))
        dist = np.array((max - origin_new)*self.n[0],dtype=int)
        new_grid = []
        for i in range(dist[0]):
                    for j in range(dist[1]):
                        for k in range(dist[2]):
                            # Coordinate of the voxel corner
                            coord = origin_new + np.array([i/(self.n[0]),j/(self.n[1]),k/(self.n[2])])
                            weight = self.cube[int(coord[0]%1)*self.n[0], int(coord[1] %1)*self.n[1], int(coord[2]%1)*self.n[2] ]
                            new_grid.append([coord[0], coord[1], coord[2], weight])
        return new_grid
    def Expand_frac(self,min,max=[]):
        if max == []:
            max = self.max_frac
        else:
            max = np.array(max)
        origin_new = self.Nearest_Coord(np.array(min))
        dist = np.array((max - origin_new)*self.n[0],dtype=int)
        new_grid = []
        for i in range(dist[0]):
                    for j in range(dist[1]):
                        for k in range(dist[2]):
                            # Coordinate of the voxel corner
                            coord = origin_new + np.array([i/(self.n[0]),j/(self.n[1]),k/(self.n[2])])
                            weight = self.cube[int(coord[0]%1)*self.n[0], int(coord[1] %1)*self.n[1], int(coord[2]%1)*self.n[2] ]
                            new_grid.append([coord[0], coord[1], coord[2], weight])
        return new_grid
    

    def around_point(self,coord):
        return
        
c = CubeFile("defc_afm1_up.cube")
grid = ProcessGrid(c)
grid.Expand(min=[-0.1,0,0])
pass

