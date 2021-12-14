from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from scipy.io import mmread
from sys import argv

mat = mmread(argv[1])
im = plt.imshow(mat.todense(), norm=LogNorm(vmin=1e-6, vmax=10))
plt.colorbar(im)
plt.savefig(argv[2])
