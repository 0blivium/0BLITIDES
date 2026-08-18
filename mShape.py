# Module: Matrix Shape Checker

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

def matrix_shape(mat,N,m1=7,m2=7,save=True):
    matrix_initialize=lambda x,m, : np.zeros((x,m),dtype=float)
    x=6*N+2
    m=m1+m2+1
    color_matrix=matrix_initialize(x,m)
    for i in range(x):
        for j in range(m):
            if (mat[i,j]!= 0):
                color_matrix[i,j]=1
    cmap = colors.ListedColormap(["white","red"])
    bounds = [0.,1.,2.]
    norm = colors.BoundaryNorm(bounds,cmap.N)
    fig,ax = plt.subplots(figsize = (10,10))   
    ax.imshow(color_matrix,cmap=cmap,norm=norm)
    ax.grid(linestyle='-',color="grey")
    ax.set_xticks(np.arange(-0.5,m,1))
    ax.set_yticks(np.arange(-0.5,x,1))
    if(save): plt.savefig("run/banded_matrix.pdf",dpi=50)
    plt.show()

    return fig

matrix=np.loadtxt("run/banded_matrix.txt")
matrix_shape(matrix,6)

####################################################################
# end of the module!
