from math import acos, pi, sqrt
from decimal import Decimal, getcontext

getcontext().prec = 30

class Vector(object):
    
    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize the zero vector'
    NO_UNIQUE_PARALLEL_COMP_MSG = 'There''s no unique vector, bra.'
    
    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Decimal(a) for a in coordinates])
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')
    
    def __add__(self, v):
        new_coordinates = [x+y for x,y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)
    
    def __sub__(self, v):
        new_coordinates = [x-y for x,y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)
        
    def __mul__(self, c):
        if isinstance(c, Vector):
            result = [x*y for x, y in zip(self.coordinates, c.coordinates)]
            return sum(result)
        else:
            return Vector([x*c for x in self.coordinates])
        
    def magnitude(self):
        coordinates_squared = [x**2 for x in self.coordinates]
        return Decimal(sqrt(sum(coordinates_squared)))
            
    def is_zero(self, tolerance=1e-10):
        return abs(self.magnitude()) < tolerance
    
    def normalized(self):
        try:
            magnitude = Decimal(self.magnitude())
            return self*(Decimal('1.0')/magnitude)
        
        except ZeroDivisionError:
            raise Exception(self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG)
    
    def angle(self, v, in_degrees=False):
        try:
            ratio = round((self.normalized())*v.normalized(),5)
            angrad = acos(ratio)
            if in_degrees:
                return angrad*Decimal(180)/pi
            else:
                return angrad
                
        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception('Cannot compute an angle with the zero vector')
            else:
                raise e
    
    def projectOn(self, b):
        try:
            normalization = b.normalized()
            return normalization*(self*normalization)
    
        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMP_MSG)
            else:
                raise e
        
    def orthogonalTo(self, b):
        try:
            projection = self.projectOn(b)
            return self.minus(projection)
        
        except Exception as e:
            if self.isOrthogonal(b) != True:
                raise Exception('Not orthogonal vectors')
            else:
                raise e
    
    def crossProductOf(self, w):
        try:
            a, b = self.coordinates, w.coordinates
            new_coordinate = [a[1]*b[2] - a[2]*b[1],
                              a[2]*b[0] - a[0]*b[2],
                              a[0]*b[1] - a[1]*b[0]]
            return Vector(new_coordinate)
            
        except ValueError as e:
            msg = str(e)
            if msg == 'need more than 2 values to unpack':
                selfInThree = Vector(self.coordinates + ('0',))
                wInThree = Vector(w.coordinates + ('0',))
                return selfInThree.crossProductOf(wInThree)
            elif (msg == 'dimension above three' or msg == 'dimension 1'):
                raise Exception(self.ONLY_DEFINE_IN_TWO_THREE_MSG)
            else:
                raise e
                
    def paraArea(self, w):
        return self.crossProductOf(w).magnitude()
        
    def triArea(self, w):
        return Decimal(0.5)*self.paraArea(w)
        
    def isParallel(self, v):
        return(self.is_zero() or v.is_zero() or self.angle(v) == 0 or self.angle(v) == pi)
    
    def isOrthogonal(self, v, tolerance=1e-10):
        return abs(self*v)<tolerance
    
    def __getitem__(self,index):
        return self.coordinates[index]
    
    def __iter__(self):
        return iter(self.coordinates)
        
    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)

    def __eq__(self, v):
        return self.coordinates == v.coordinates