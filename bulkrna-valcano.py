#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
from adjustText import adjust_text

    
    df_sorted = df.sort_values(by='padj', ascending=False)
    plt.scatter(x=df['log2FoldChange'],y=df['padj'].apply(lambda x:-np.log10(x)),s=100,label="Not significant",color = 'grey',alpha=0.35,edgecolor='none')

    down_df = df[(df['log2FoldChange']<=-1)&(df['padj']<=0.05)]
    up_df = df[(df['log2FoldChange']>=1)&(df['padj']<=0.05)]

    plt.scatter(x=down_df['log2FoldChange'],y=down_df['padj'].apply(lambda x:-np.log10(x)),s=100,label="Down-regulated",color="teal",alpha=0.55,edgecolor='none')
    plt.scatter(x=up_df['log2FoldChange'],y=up_df['padj'].apply(lambda x:-np.log10(x)),s=100,label="Up-regulated",color="salmon",alpha=0.55,edgecolor='none')
    
    
    down = df[(df['log2FoldChange']<=-1)&(df['padj']<=0.05)][:10]
    up = df[(df['log2FoldChange']>=1)&(df['padj']<=0.05)][:10]
    
    texts_up=[]
    for i,r in up.iterrows():
        if r['padj']== 0:
            texts_up.append(plt.text(x=r['log2FoldChange'],y=-np.log10(1e-308),s=i,fontsize=12))
        else:
            texts_up.append(plt.text(x=r['log2FoldChange'],y=-np.log10(r['padj']),s=i,fontsize=12))
    
    texts_down=[]
    for i,r in down.iterrows():
        if r['padj']==0:
            texts_down.append(plt.text(x=r['log2FoldChange'],y=-np.log10(1e-307),s=i,fontsize=12))
        else:
            texts_down.append(plt.text(x=r['log2FoldChange'],y=-np.log10(r['padj']),s=i,fontsize=12))
    
    adjust_text(texts_up+texts_down,arrowprops=dict(arrowstyle="-", color='black', lw=0.5))


    plt.figure(12,8)
    plt.xlabel("logFC")
    plt.ylabel("-logFDR")
    plt.ylim(-50,350)
    plt.axvline(-1,color="grey",linestyle="--")
    plt.axvline(1,color="grey",linestyle="--")
    plt.axhline(1,color="grey",linestyle="--")
    plt.legend(loc='upper left')

