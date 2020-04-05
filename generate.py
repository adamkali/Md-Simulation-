from sys import argv
from random import uniform

__file__ = open("input.xy", "w")

box = float(argv[1])

for i in range(50):
    x, y, z = uniform(-box, box), uniform(-box,box), uniform(-box,box)
    print(x,y,z, file = __file__)

