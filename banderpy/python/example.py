from bandersnatch import Point, Scalar
import copy
import time

# constructor
s = Scalar()
s.from_u64(0x123)
print(s)

# get the generator
g = Point(generator = True)
print("generator:", g)

# get a random point
p = Point()
print("random point:", p)
p1 = copy.deepcopy(p)
p2 = copy.deepcopy(p)

# get a random scalar
s = Scalar()
print("random scalar:", s)

# scalar multiplication via double-then-add
t1 = time.process_time()
print(t1)
for i in range(1000):
    p1.mul(s)
t1 = time.process_time() - t1
print("double-then-add time:", t1)

# scalar multiplication via glv
t2 = time.process_time()
for i in range(1000):
    p2.glv(s)
t2 = time.process_time() - t2
print("glv time:", t2)

assert p1 == p2

# multi-base-multiplication
r = Point()
bases = [Point() for _ in range (1000)]
scalars = [Scalar() for _ in range (1000)]
t3 = time.process_time()
Point.msm(r, bases, scalars)
t3 = time.process_time() - t3
print("msm time:", t3)

# doubling
p = Point()
p1 = copy.deepcopy(p)
p2 = copy.deepcopy(p)
p1.double()
p2.add(p)

assert p1 == p2

# non-compressed form
print(p.p)
# compressed form
print(p.serialize())

# serdes a point
b = p.serialize()
p2 = Point()
p2.deserialize(b)

assert p == p2

# serdes a scalar
s = Scalar()
b = s.serialize()
s2 = Scalar()
s2.deserialize(b)

assert s == s2