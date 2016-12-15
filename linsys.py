from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane
from hyperplane import Hyperplane

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)
            
    def do_gaussian_elimination_and_parametrization(self):
        rref = self.compute_rref()
        rref.raise_excepion_if_contradictory_equation()

        direction_vectors = rref.extract_direction_vectors_for_parametrization()  # NOQA
        basepoint = rref.extract_basepoint_for_parametrization()

        return Parametrization(basepoint, direction_vectors)

    def extract_direction_vectors_for_parametrization(self):
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        free_variable_indices = set(range(num_variables)) - set(pivot_indices)

        direction_vectors = []

        for free_var in free_variable_indices:
            vector_coords = [0] * num_variables
            vector_coords[free_var] = 1
            for index, plane in enumerate(self.planes):
                pivot_var = pivot_indices[index]
                if pivot_var < 0:
                    break
                vector_coords[pivot_var] = -plane.normal_vector[free_var]

            direction_vectors.append(Vector(vector_coords))

        return direction_vectors

    def extract_basepoint_for_parametrization(self):
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()

        basepoint_coords = [0] * num_variables

        for index, plane in enumerate(self.planes):
            pivot_var = pivot_indices[index]
            if pivot_var < 0:
                break
            basepoint_coords[pivot_var] = plane.constant_term

        return Vector(basepoint_coords)

    def raise_excepion_if_contradictory_equation(self):
        for plane in self.planes:
            try:
                plane.first_nonzero_index(plane.normal_vector)

            except Exception as e:
                if str(e) == 'No nonzero elements found':
                    constant_term = MyDecimal(plane.constant_term)
                    if not constant_term.is_near_zero():
                        raise Exception(self.NO_SOLUTIONS_MSG)

                else:
                    raise e

    def raise_excepion_if_too_few_pivots(self):
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        num_pivots = sum([1 if index >= 0 else 0 for index in pivot_indices])
        num_variables = self.dimension

        if num_pivots < num_variables:
            raise Exception(self.INF_SOLUTIONS_MSG)

    def compute_solution(self):
        """
        Return the solutions of the System
        # if there is 0 = k there are no solutions
        # if there are less planes than dimensions there are infinite solutions
        # if there are more than 1 pivot variable there are infinite solutions
        # otherwise there is a single solution
        """
        try:
            return self.do_gaussian_elimination_and_parametrization()

        except Exception as e:
            if (str(e) == self.NO_SOLUTIONS_MSG):
                return str(e)
            else:
                raise e
                
#    def system_solutions(self):
#        
#        rref = self.compute_rref()
#        
#        print('\nRREF spits out {} \n'.format(rref))
#        indices = rref.indices_of_first_nonzero_terms_in_each_row()
#        if (-1 in indices):
#            for i, j in enumerate(indices):
#                constant_is_zero = MyDecimal(rref.planes[i].constant_term).is_near_zero()
#                if j == -1 and not constant_is_zero:
#                    print ("No solution")
#                    break
                # Currently flawed must add parameterization
#                elif i == len(indices) -1:
#                    print ("Infinte solutions")
#        else:
#            print ("The solution is : ")
#            print (rref)
        
    def compute_rref(self):
        tf = self.compute_triangular_form()
        num_equations =len(tf)
        try:
            tf[num_equations-1].normal_vector[num_equations-1]
        except IndexError:
            num_equations -= 1
        for num_row in range(num_equations-1,-1, -1):
            leading_row_is_zero = MyDecimal(tf[num_row].normal_vector[num_row]).is_near_zero()
            if tf[num_row].normal_vector[num_row] != 1 and not leading_row_is_zero:
                    tf[num_row] *= Decimal(1)/tf[num_row].normal_vector[num_row]
            # Loop that eliminates nonzeroes below current row.
            for index in range(num_row-1, -1,-1):
                bottom_rows_check = MyDecimal(tf[index].normal_vector[num_row]).is_near_zero()
                if not bottom_rows_check and not leading_row_is_zero:
                    # /tf[num_row].normal_vector[num_row]
                    coeff = tf[index].normal_vector[num_row]
                    tf.add_multiple_times_row_to_row(-coeff, num_row, index)
        return tf
    
    def compute_triangular_form(self):
        triform = deepcopy(self)
        num_equations = len(triform)
        for num_row in range(num_equations):
            # Loop that swaps rows for leading coeffcients if possible.
            for index in range(num_row + 1, num_equations):
                if triform[num_row].normal_vector[num_row] == 0 and triform[index].normal_vector[num_row] != 0:
                    triform.swap_rows(index-1, index)
                    break
            # Loop that eliminates nonzeroes below current row.
            for index in range(num_row + 1, num_equations):
                bottom_rows_check = MyDecimal(triform[index].normal_vector[num_row]).is_near_zero()
                if not  bottom_rows_check:
                    coeff = triform[index].normal_vector[num_row]/triform[num_row].normal_vector[num_row]
                    triform.add_multiple_times_row_to_row(-coeff, num_row, index)
        return triform

    def swap_rows(self, row1, row2):
        self[row1], self[row2] = self[row2], self[row1]


    def multiply_coefficient_and_row(self, coefficient, row):
        #redundant but added to pass course test cases
        self[row] *= coefficient


    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        plane_add, plane_added = self[row_to_add], self[row_to_be_added_to]
        self[row_to_be_added_to] = plane_add * coefficient + plane_added

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices


    def __len__(self):
        return len(self.planes)


    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        ret = 'Linear system:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps

class Parametrization(object):

    BASEPT_AND_DIR_VECTORS_MUST_BE_IN_SAME_DIM = (
        'The basepoint and direction vectors should all live in the same '
        'dimension')

    def __init__(self, basepoint, direction_vectors):

        self.basepoint = basepoint
        self.direction_vectors = direction_vectors
        self.dimension = self.basepoint.dimension

        try:
            for v in direction_vectors:
                assert v.dimension == self.dimension

        except AssertionError:
            raise Exception(self.BASEPT_AND_DIR_VECTORS_MUST_BE_IN_SAME_DIM)

    def __str__(self):

        output = ''
        for coord in range(self.dimension):
            output += 'x_{} = {} '.format(coord + 1,
                                          round(self.basepoint[coord], 3))
            for free_var, vector in enumerate(self.direction_vectors):
                output += '+ {} t_{}'.format(round(vector[coord], 3),
                                             free_var + 1)
            output += '\n'
        return output

p1 = Plane(Vector([0.786, 0.786, 0.588]), -0.714)
p2 = Plane(Vector([-0.131, -0.131, 0.244]), 0.319)

system = LinearSystem([p1, p2])
print (system.compute_solution())


p1 = Plane(Vector([8.631, 5.112, -1.816]), -5.113)
p2 = Plane(Vector([4.315, 11.132, -5.27]), -6.775)
p3 = Plane(Vector([-2.158, 3.01, -1.727]), -0.831)

system = LinearSystem([p1, p2, p3])
print (system.compute_solution())

p1 = Plane(Vector([0.935, 1.76, -9.365]), -9.955)
p2 = Plane(Vector([0.187, 0.352, -1.873]), -1.991)
p3 = Plane(Vector([0.374, 0.704, -3.746]), -3.982)
p4 = Plane(Vector([-0.561, -1.056, 5.619]), 5.973)

system = LinearSystem([p1, p2, p3, p4])
print (system.compute_solution())


