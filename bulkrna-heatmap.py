#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import seaborn as sns

    results = []
    for gene in gene_list:
        indices = df[df['symbol'] == gene].index.tolist()
        for index in indices:
            results.append({'Symbol': gene, 'Index': index})
    GOI_df = pd.DataFrame(results)
    GOI_df = GOI_df.set_index('Index')
    
    dds_GOI = dds_me[:,GOI_df.index]
    heatmap_data = pd.DataFrame(dds_GOI.layers['log1p'].T,
                                index=dds_GOI.var_names, columns=dds_GOI.obs_names)

    heatmap_data = heatmap_data[order_req]
    heatmap_data['symbol'] = df['symbol']
    heatmap_data.set_index('symbol', inplace=True)

    hm = sns.clustermap(heatmap_data, cmap="coolwarm", col_cluster=False, row_cluster= False,linewidths=0.5, linecolor='white',yticklabels=True,mask=False)

    hm.savefig(f"{file_name}.png", dpi=300)
    
    plt.show()

