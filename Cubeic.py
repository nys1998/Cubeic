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

    
    def get_coordinates(self):
        """Generate physical coordinates for all grid points."""
        grid_coords = []
        shape = self.n
        axis = self.axis
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    # Coordinate of the voxel corner
                    coord = self.origin + (i * axis[0]/shape[0] +
                                            j * axis[1]/shape[1] + k * axis[2]/shape[2])
                    weight = self.data[i, j, k]
                    grid_coords.append([coord[0], coord[1], coord[2], weight])
        return grid_coords

class ProcessGrid():
    '''
    Process the Cube grid in xyzw format
    '''
    def __init__(self, cube:CubeFile):
        self.data = cube.get_coordinates()
        self.axis = cube.axis
        self.n = cube.n

    def Expand(self,x,y,z):
        pass
        
c = CubeFile("defc_afm1_up.cube")
xyz = c.get_coordinates()
pass

