import numpy as np
import matplotlib.pyplot as plt
import _astar

grid = np.array([[2,1,0,0,0],
                 [0,5,6,0,0],
                 [0,1,0,2,0],
                 [0,0,1,0,0]], dtype=np.float64)
grid2 = np.array([[0,0,0,0,0],
                  [0,0,0,0,0],
                  [0,0,0,0,0],
                  [0,0,0,0,0]], dtype=np.int)

grid3 = np.zeros((1200,1200))
cost = _astar.astar(1,1,3,3,grid2)
#cost = _astar.astar(100,100,900,900,grid3)
print cost

plt.imshow(cost, interpolation='None')
plt.show()
