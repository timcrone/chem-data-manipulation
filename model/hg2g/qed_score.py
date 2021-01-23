import sys
sys.path.append('../')
from props import *

for line in sys.stdin:
    x,y = line.split()
    if y == "None": y = None
    sim2D = similarity(x, y)
    try:
        print(x, y, sim2D,qed(x), qed(y))
    except Exception as e:
        print(x, y, sim2D,0.0, 0.0)
