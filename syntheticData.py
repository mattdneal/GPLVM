import numpy as np
data_temp = np.genfromtxt("/media/matt/Data1/Dropbox/R_packages/GPLVM/syntheticData.csv")
data_temp = data_temp[1:]
data = np.reshape(data_temp, (1000,5,5), "F")
#data_flat = np.reshape(data_temp, (1000,25), "F")
subset_size = 200
indices = np.random.choice(range(1000), subset_size)
data = data[indices, :, :]
data_flat = np.reshape(data, (subset_size, -1))

import GPflow as gp

pca = gp.gplvm.PCA_reduce(data_flat, 2)

%matplotlib inline
import matplotlib.pyplot as plt
def pltScatter(data):
    plt.scatter(data[:, 0], data[:, 1])
#
pltScatter(pca)

#gplvm = gp.gplvm.GPLVM(data_flat, 2)
#gplvm.optimize(disp=True, maxiter=100)

sgplvm = gp.gplvm.SGPLVM(data, 2)

from copy import deepcopy
X_latent_init = deepcopy(sgplvm.X_latent.value)
pltScatter(X_latent_init)
sgplvm.optimize(maxiter=100, disp=True)
X_latent_final = deepcopy(sgplvm.X_latent.value)
pltScatter(X_latent_final)
