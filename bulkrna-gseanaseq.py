#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from gseapy import GSEA

    gs = GSEA(data=df,
              gene_sets="mh.all.v2023.2.Mm.symbols.gmt",
              classes=cls_path,
              permutation_type='phenotype',
              permutation_num=1000,
              outdir=None,
              method='signal_to_noise',
              threads=4,
              seed=8)
    gs.pheno_pos = pheno_pos
    gs.pheno_neg = pheno_neg
    gs.run()


    out = []
    for term in list(gs.results):
        if gs.results[term]['fdr'] < 0.05:
            out.append([
                term,
                gs.results[term]['fdr'],
                gs.results[term]['es'],
                gs.results[term]['nes']
            ])

    out_df = pd.DataFrame(out, columns=['Term', 'fdr', 'es', 'nes'])
    out_df = out_df.sort_values('fdr').reset_index(drop=True)


# In[ ]:


from gseapy import dotplot

ax_1 = dotplot(gs.res2d,
             column="FDR q-val",
             ofname = None,
             title = '####',
             cmap=plt.cm.viridis,
             size=4,legend_size = 1,
             figsize=(4,5), cutoff=1)

