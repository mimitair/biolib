###############
### CLASSES ###
###############
from dataclasses import dataclass

@dataclass
class Point:
    """
    This class represents a point in 3D space.
    """
    x: float
    y: float
    z: float
    
    def distance(self, other: 'Point'):
        """
        Calculate Euclidian distance between this Point and another Point object.
        """
        # DEFENSIVE CHECK:
        if not isinstance(other, Point):
            raise TypeError("'other' must also be a Point object to calculate distance")
        
        # Import square root:
        from math import sqrt  

        # Calculate the distances between each coordinate in its respective axis (x, y and z)
        x_dist: float = (self.x - other.x)**2
        y_dist: float = (self.y - other.y)**2
        z_dist: float = (self.z - other.z)**2
        
        return sqrt(x_dist + y_dist + z_dist)



#################
### FUNCTIONS ###
#################
def printTriad(res1: str, 
               res2: str, 
               res3: str, 
               atom1: str,
               atom2: str,
               atom3: str,
               min_dist_1_2: float, 
               min_dist_2_3: float,
               min_dist_3_1: float,
               max_dist_1_2: float,
               max_dist_2_3: float,
               max_dist_3_1: float):
    
    """
    Function to print a given amino acid triad in a pretty format :)
    """
    print(f"              {res1}({atom1}) ")
    print("               /  |")
    print("              /   |")
    print("             /    |")
    print("            /     |")
    print(f"   ({min_dist_1_2} - {max_dist_1_2}) /      | ({min_dist_3_1} - {max_dist_3_1})")
    print("          /       |")
    print("         /        |")
    print("        /         |")
    print("       /          |")
    print(f"{res2}({atom2}) ------- {res3}({atom3})")
    print(f"         ({min_dist_2_3} - {max_dist_2_3})")
 
