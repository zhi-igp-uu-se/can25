#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sklearn
import scanpy as sc
import pandas as pd
from sklearn import datasets, preprocessing, decomposition, model_selection
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


scaler = StandardScaler()
scaled_tpm = scaler.fit_transform(tpm)
pca = PCA(n_components=2) 
pca_result = pca.fit_transform(scaled_tpm)

pca_tpm = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
pca_tpm_index=pca_tpm

pca_tpm_index["Cell line"]=metadata["Cell line"]
pca_tpm_index["Batch"]=metadata["Batch"]
pca_tpm_index["Batch"]= pca_tpm_index["Batch"].astype(str)

color_map = {'CCR7 PDPN': 'red', 'CCR7': 'blue', 'X2 PDPN': 'green','X2':"purple"}
marker_map = {'1': 'o', '2': 's'}

plt.figure(figsize=(8, 6))
texts = []

for i in range(len(pca_tpm)):
    genotype = pca_tpm_index.iloc[i]['Cell line']
    batch = pca_tpm_index.iloc[i]['Batch']
    
    plt.scatter(pca_result[i, 0], pca_result[i, 1], color=color_map[genotype], marker=marker_map[batch])

legend_elements = [Line2D([0], [0], marker=marker_map[batch], color='w', label='Batch: ' + batch,
                          markerfacecolor='black', markersize=10) for batch in marker_map]
legend_elements += [Line2D([0], [0], color=color_map[genotype], lw=5, label='Cell line: ' + genotype) for genotype in color_map]

plt.legend(handles=legend_elements)
plt.title('PCA of TPM')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.show()

