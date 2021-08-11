from libbanderpy import *
from copy import deepcopy

class Point():

    def __init__(self, generator = False):
        if generator==False:
            self.p = random_point_rust()
        else: 
            self.p = get_generator_rust()
        
    def __str__(self):
        return point_to_string_rust(self.p)

    def __eq__(self, other):
        return self.p == other.p
    
    def add(self, other):
        self.p = add_rust(self.p, other.p)
        return self

    def double(self):
        self.p = double_rust(self.p)
        return self

    def mul(self, scalar):
        if isinstance(scalar, int):
            self.p = mul_rust(self.p, Scalar().from_int(scalar).s)
        else: 
            self.p = mul_rust(self.p, scalar.s)
        return self

    def glv(self, scalar):
        if isinstance(scalar, int):
            self.p = glv_rust(self.p, Scalar().from_int(scalar).s)
        else: 
            self.p = glv_rust(self.p, scalar.s)
        return self

    def msm(self, points, scalars):
        self.p = msm_rust([x.p for x in points],
                          [Scalar().from_int(x).s if isinstance(x, int) else x.s for x in scalars])
        return self

    def dup(self):
        return deepcopy(self)

    def serialize(self):
        return bytes(point_serialize_rust(self.p))

    def deserialize(self, b):
        s = [e for e in b]
        self.p = point_deserialize_rust(s)
        return self


class Scalar():
    def __init__(self):
        self.s = random_scalar_rust()

    def dup(self):
        return deepcopy(self)

    def __str__(self):
        return scalar_to_string_rust(self.s)

    def __eq__(self, other):
        return self.s == other.s

    def serialize(self):
        return bytes(self.s)

    def deserialize(self, b):
        s = [e for e in b]
        self.s = s
        return self

    def from_int(self, b):
        self.s = [e for e in b.to_bytes(32, "little")]
        return self