import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

n_ct_bins = 80  # adjust as needed

n_cths = np.arange(1, n_ct_bins+1) # The bin choices

df = pd.read_csv(f'chi2_all64_nctbins{n_ct_bins}.csv')

ith_layer = 64  # Isa change this to explore different layers, each ith_layerums correspond to a perturbed layer 

if ith_layer > 64:
    print("Error: Maximum number of layers is 64")

else:

    col = ith_layer - 1
    plt.figure(figsize=(8, 5))
    plt.plot(n_cths, df.iloc[:, col], marker='o', label=f' layer {ith_layer}')
    plt.xlabel('ct bins')
    plt.ylabel(rf'$\chi^2$')
    plt.title(rf'$\chi^2$ vs nbins — ith layer {ith_layer}')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()