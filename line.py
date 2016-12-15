from decimal import Decimal, getcontext
from vector import Vector

getcontext().prec = 30


class Line(object):

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 2

        if not normal_vector:
            all_zeros = ['0']*self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = normal_vector

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        self.set_basepoint()


    def set_basepoint(self):
        try:
            n = self.normal_vector
            c = self.constant_term
            basepoint_coords = ['0']*self.dimension

            initial_index = Line.first_nonzero_index(n)
            initial_coefficient = Decimal(n[initial_index])

            basepoint_coords[initial_index] = c/initial_coefficient
            self.basepoint = Vector(basepoint_coords)
            if not isinstance(self.normal_vector, Vector):
                self.normal_vector = Vector(self.normal_vector)

        except Exception as e:
            if str(e) == Line.NO_NONZERO_ELTS_FOUND_MSG:
                self.basepoint = None
            else:
                raise e


    def __str__(self):

        num_decimal_places = 3

        def write_coefficient(coefficient, is_initial_term=False):
            coefficient = round(coefficient, num_decimal_places)
            if coefficient % 1 == 0:
                coefficient = int(coefficient)

            output = ''

            if coefficient < 0:
                output += '-'
            if coefficient > 0 and not is_initial_term:
                output += '+'

            if not is_initial_term:
                output += ' '

            if abs(coefficient) != 1:
                output += '{}'.format(abs(coefficient))

            return output

        n = self.normal_vector

        try:
            initial_index = Line.first_nonzero_index(n)
            terms = [write_coefficient(n[i], is_initial_term=(i==initial_index)) + 'x_{}'.format(i+1)
                     for i in range(self.dimension) if round(n[i], num_decimal_places) != 0]
            output = ' '.join(terms)

        except Exception as e:
            if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
                output = '0'
            else:
                raise e

        constant = round(self.constant_term, num_decimal_places)
        if constant % 1 == 0:
            constant = int(constant)
        output += ' = {}'.format(constant)

        return output


    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Line.NO_NONZERO_ELTS_FOUND_MSG)
    
    def isParallelTo(self, v):
        # Take the base points, form a vector, and check to see if
        # the normal of the new vector is the same as the normal of sample vectors
        return self.normal_vector.isParallel(v.normal_vector)
    
    def __eq__(self, v):
        if self.normal_vector.is_zero():
            if not v.normal_vector.is_zero():
                return False
            else:
                diff = self.constant_term - v.constant_term
                return MyDecimal(diff).is_near_zero()
        elif v.normal_vector.is_zero():
            return False
        
        # Check if parallel
        if not self.isParallelTo(v):
            return False
        
        l1Base = self.basepoint
        l2Base = v.basepoint
        # Create a new line between a point on each
        temp_vector = l1Base-l2Base
        # If it's the same line, then the new line will be orthogonal 
        # to both lines
        return self.normal_vector.isOrthogonal(temp_vector) 
        return False
    
    def intersectionWith(self, v):
        try:
            if self.isParallelTo(v):
                return self
            # Intialization of line 1 and Line 2
            l1, l2 = self.normal_vector.coordinates, v.normal_vector.coordinates
            # Using the derived formula x = (A*k2-C*k1)/(AD-BC) and y = (D*k1-B*k2)/(AD-BC)
            v_normalized = Vector([l2[1], -l2[0]])
            AD_BC = self.normal_vector*v_normalized
            A, B = l1
            C, D = l2
            k1, k2 = self.constant_term, v.constant_term
            y, x = A*k2-C*k1, D*k1-B*k2
            return Vector([x,y])*(Decimal(1)/AD_BC)
                
        except ZeroDivisionError:
            print(self == v)
            if self == v:
                return self
            else:
                return None

class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps