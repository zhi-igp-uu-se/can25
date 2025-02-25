#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import scanpy as sc
import pandas as pd


# In[ ]:


# umap


# In[ ]:


sc.pp.neighbors(adata, n_pcs= 50)
sc.tl.umap(adata, n_components = 2)


# In[ ]:


subtypes = adata.obs['subtype'].unique()
colors = sc.pl.palettes.vega_10[:len(subtypes)] 
subtype_color_map = dict(zip(subtypes, colors))


# In[ ]:


umap_df = pd.DataFrame(adata_cell.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata_cell.obs.index)
umap_df['subtype'] = adata_cell.obs['subtype']
umap_df['cell_line'] = adata_cell.obs['sample']


# In[ ]:


fig, ax = plt.subplots(figsize=(10, 8))
for subtype, color in subtype_color_map.items():
    subset = umap_df[umap_df['subtype'] == subtype]
    ax.scatter(subset['UMAP1'], subset['UMAP2'], s=5, label=subtype, c=color, alpha=0.7)

for cell_line, (x, y) in cell_line_centers.iterrows():
    ax.text(x, y, cell_line, fontsize=10, fontweight='bold', ha='center', color='black')


ax.legend(
    title="Subtype",
    loc="lower right",       
    bbox_to_anchor=(1.25, 0.1),  
    frameon=True            
)

ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
ax.set_xticks([])
ax.set_yticks([])
ax.grid(False)
plt.show()


# In[ ]:


sc.pl.umap(
    adata,
    color=['PDPN'],          
    cmap='Reds', 
    layer='log1p_norm',      
    vmin=0                  
)
figsize=(10, 8)

