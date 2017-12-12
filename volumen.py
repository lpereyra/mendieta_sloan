import numpy as np
from scipy.spatial import ConvexHull

bintype = np.dtype([
        ("x", np.float32),
        ("y", np.float32),
        ("z", np.float32),
                  ])

binario = open("sloan_0.10_xyz.bin","rb");
N = np.fromfile(binario, dtype=np.int32, count=1); N = N[0]
print N
data = np.fromfile(binario, dtype=bintype)

x = data["x"]
y = data["y"]
z = data["z"]

points = zip(x,y,z)
hull = ConvexHull(points)

print "volumen ConvexHull: ",hull.volume
