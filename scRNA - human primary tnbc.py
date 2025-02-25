#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import scanpy as sc


# In[ ]:


sc.tl.leiden(adata, key_added="leiden_res0_1", resolution=0.1)


# In[ ]:


# umap
sc.pl.umap(
    adata,
    color=["leiden_res0_1"],
    legend_loc="on data",
)


# In[ ]:


# dot plot
sc.tl.dendrogram(adata, groupby='leiden_res0_1')
sc.pl.dotplot(adata, markers, groupby='leiden_res0_1', dendrogram=True, layer = 'log1p_norm',swap_axes=True,var_group_rotation
=90)


# In[ ]:


# gsea
data_dense = adata.layers['raw_counts'].todense()
data = pd.DataFrame(data_dense, 
                    index=adata.obs_names, 
                    columns=adata.var_names)
 
data.columns = data.columns.astype(str).str.upper()
data.index = data.index.astype(str)

data = data.T


# In[ ]:


gene_sets = gp.parser.read_gmt('h.all.v2023.2.Hs.symbols.gmt')
res_cluster = gp.gsea(
    data=data,  
    gene_sets=gene_sets,  
    cls=cls_cluster,  
    permutation_num=1000,
    permutation_type='phenotype',
    
    outdir=None,
    method='s2n',  
    min_size=5,    
    max_size=5000, 
    threads=16
)
res_cluster.pheno_pos = "0"
res_cluster.pheno_neg = "1"
res_cluster.run()


# In[ ]:


ax = barplot(selected_terms, column='FDR q-val', title='Cluster 0 vs Cluster 1', figsize=(6, 4))

