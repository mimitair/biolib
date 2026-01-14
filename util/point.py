from math import sqrt  # For square root calculation

class Point:
    """
    This class represents a point in 3D space.
    """
    
    def __init__(self, x: float, y: float, z: float):
        ### DEFENSIVE CHECKS ###
        if not isinstance(x, float) & isinstance(y, float) & isinstance(z, float) :
            raise TypeError("Point class only accepts floats as input.")
        
        ### INIT ###
        self.x: float = x
        self.y: float = y
        self.z: float = z
    
    
    def distance(self, other: 'Point'):
        """
        Calculate Euclidian distance between this Point and another Point object.
        """
        # DEFENSIVE CHECK:
        if not isinstance(other, Point):
            raise TypeError("'other' must also be a Point object to calculate distance")
        
        # Calculate the distances between each coordinate in its respective axis (x, y and z)
        x_dist: float = (self.x - other.x)**2
        y_dist: float = (self.y - other.y)**2
        z_dist: float = (self.z - other.z)**2
        
        return sqrt(x_dist + y_dist + z_dist)
