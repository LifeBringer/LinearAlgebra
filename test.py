from linsys import LinearSystem
from plane import Plane
from line import Line
from vector import Vector
from decimal import Decimal, getcontext

#Test case from Line
A, B = Line(Vector([4.046, 2.836]), 1.21), Line(Vector([10.115, 7.09]), 3.025)
C, D = Line(Vector([7.204, 3.182]), 8.68), Line(Vector([8.172, 4.114]), 9.883)
E, F = Line(Vector([1.182,5.562]), 6.744), Line(Vector([1.773, 8.343]), 9.525)

print('\n')
print(A.intersectionWith(B))
print(A == B)

print('\n')
print(C.intersectionWith(D))
print(C == D)

print('\n')
#print(E.intersectionWith(F))
print(E == F)

#Test case from Plane
A, B = Plane(Vector(['-0.412', '3.806', '0.728']), '-3.46'), Plane(Vector(['1.03', '-9.515', '-1.82']), '8.65')
C, D = Plane(Vector(['2.611', '5.528', '0.283']), '4.6'), Plane(Vector(['7.715', '8.306', '5.342']), '3.76')
E, F = Plane(Vector(['-7.926', '8.625', '-7.212']), '-7.952'), Plane(Vector(['-2.642', '2.875', '-2.404']), '2.443')

print('\n')
print(A.isParallelTo(B))
print(A == B)

print('\n')
print(C.isParallelTo(D))
print(C == D)

print('\n')
print(E.isParallelTo(F))
print(E == F)

#Linear System Operation Test Cases
p0 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p1 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
p2 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
p3 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')

s = LinearSystem([p0,p1,p2,p3])
s.swap_rows(0,1)
if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
    print('test case 1 failed')

s.swap_rows(1,3)
if not (s[0] == p1 and s[1] == p3 and s[2] == p2 and s[3] == p0):
    print( 'test case 2 failed')

s.swap_rows(3,1)
if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
    print( 'test case 3 failed')

s.multiply_coefficient_and_row(1,0)
if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
    print( 'test case 4 failed')

s.multiply_coefficient_and_row(-1,2)
if not (s[0] == p1 and
        s[1] == p0 and
        s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
        s[3] == p3):
    print( 'test case 5 failed')

s.multiply_coefficient_and_row(10,1)
if not (s[0] == p1 and
        s[1] == Plane(normal_vector=Vector(['10','10','10']), constant_term='10') and
        s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
        s[3] == p3):
    print( 'test case 6 failed')

s.add_multiple_times_row_to_row(0,0,1)
if not (s[0] == p1 and
        s[1] == Plane(normal_vector=Vector(['10','10','10']), constant_term='10') and
        s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
        s[3] == p3):
    print( 'test case 7 failed')

s.add_multiple_times_row_to_row(1,0,1)
if not (s[0] == p1 and
        s[1] == Plane(normal_vector=Vector(['10','11','10']), constant_term='12') and
        s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
        s[3] == p3):
    print( 'test case 8 failed')

s.add_multiple_times_row_to_row(-1,1,0)
if not (s[0] == Plane(normal_vector=Vector(['-10','-10','-10']), constant_term='-10') and
        s[1] == Plane(normal_vector=Vector(['10','11','10']), constant_term='12') and
        s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
        s[3] == p3):
    print( 'test case 9 failed')

#Linear System Triangular Form Test Cases
p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0','1','1']), constant_term='2')
s = LinearSystem([p1,p2])
t = s.compute_triangular_form()
if not (t[0] == p1 and
        t[1] == p2):
    print ('test case 1 failed')

p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1','1','1']), constant_term='2')
s = LinearSystem([p1,p2])
t = s.compute_triangular_form()
if not (t[0] == p1 and
        t[1] == Plane(constant_term='1')):
    print ('test case 2 failed')

p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
p4 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')
s = LinearSystem([p1,p2,p3,p4])
t = s.compute_triangular_form()
if not (t[0] == p1 and
        t[1] == p2 and
        t[2] == Plane(normal_vector=Vector(['0','0','-2']), constant_term='2') and
        t[3] == Plane()):
    print ('test case 3 failed')

p1 = Plane(normal_vector=Vector(['0','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1','-1','1']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1','2','-5']), constant_term='3')
s = LinearSystem([p1,p2,p3])
t = s.compute_triangular_form()
if not (t[0] == Plane(normal_vector=Vector(['1','-1','1']), constant_term='2') and
        t[1] == Plane(normal_vector=Vector(['0','1','1']), constant_term='1') and
        t[2] == Plane(normal_vector=Vector(['0','0','-9']), constant_term='-2')):
    print ('test case 4 failed')
    
# Test for RREF
print('\n Test for RREF:')

p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0','1','1']), constant_term='2')
s = LinearSystem([p1,p2])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term='-1') and
        r[1] == p2):
    print ('test case 1 failed')

p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1','1','1']), constant_term='2')
s = LinearSystem([p1,p2])
r = s.compute_rref()
if not (r[0] == p1 and
        r[1] == Plane(constant_term='1')):
    print ('test case 2 failed')

p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
p4 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')
s = LinearSystem([p1,p2,p3,p4])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term='0') and
        r[1] == p2 and
        r[2] == Plane(normal_vector=Vector(['0','0','-2']), constant_term='2') and
        r[3] == Plane()):
    print ('test case 3 failed')

p1 = Plane(normal_vector=Vector(['0','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1','-1','1']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1','2','-5']), constant_term='3')
s = LinearSystem([p1,p2,p3])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term=Decimal('23')/Decimal('9')) and
        r[1] == Plane(normal_vector=Vector(['0','1','0']), constant_term=Decimal('7')/Decimal('9')) and
        r[2] == Plane(normal_vector=Vector(['0','0','1']), constant_term=Decimal('2')/Decimal('9'))):
    print ('test case 4 failed')
    
# Test for solution
plane_1 = Plane(normal_vector=Vector([5.862, 1.178, -10.366]), constant_term='-8.15')
plane_2 = Plane(normal_vector=Vector([-2.931, -0.589, 5.183]), constant_term='-4.075')

print('System 1')
lin_sys_1 = LinearSystem([plane_1, plane_2])
print(lin_sys_1.compute_solution())


plane_3 = Plane(normal_vector=Vector([8.631, 5.112, -1.816]), constant_term='-5.113')
plane_4 = Plane(normal_vector=Vector([4.315, 11.132, -5.27]), constant_term='-6.775')
plane_5 = Plane(normal_vector=Vector([-2.158, 3.01, -1.727]), constant_term='-0.831')

print('System 2')
lin_sys_2 = LinearSystem([plane_3, plane_4, plane_5])
print(lin_sys_2.compute_solution())


plane_6 = Plane(normal_vector=Vector([5.262, 2.739, -9.878]), constant_term='-3.441')
plane_7 = Plane(normal_vector=Vector([5.111, 6.358, 7.638]), constant_term='-2.152')
plane_8 = Plane(normal_vector=Vector([2.016, -9.924, -1.367]), constant_term='-9.278')
plane_9 = Plane(normal_vector=Vector([2.167, -13.543, -18.883]), constant_term='-10.567')

print('System 3')
lin_sys_3 = LinearSystem([plane_6, plane_7, plane_8, plane_9])
print(lin_sys_3.compute_solution())