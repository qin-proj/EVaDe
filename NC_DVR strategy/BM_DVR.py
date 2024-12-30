import pyreadr
import pandas as pd
import scipy
from scipy import stats
import numpy as np
from matplotlib import pyplot
exp=pd.read_csv(f'/data1/qintian/science_repeat/data_mouse/GSE214390/bm_data/mean_expression_nmr.csv',index_col=0)
exp2=pd.read_csv(f'/data1/qintian/science_repeat/data_mouse/GSE214390/bm_data/mean_expression_mm2.csv',index_col=0)

cp_list=["immature.B.cell","erythroid","metamyelocyte","promyelocyte","pre.T.cell","myelocyte","B.cell","RBC","monocyte","neutrophil"]
cp_list2=["erythroid","metamyelocyte","promyelocyte","pre.T.cell","myelocyte","B.cell","RBC","monocyte","neutrophil"]
df0 = pd.read_csv("/data1/qintian/science_repeat/data_mouse/GSE214390/result/sum_frame/immature.B.cell.csv",header=0)
#0:sp 1:sam 2:Residuals
df0=df0.rename(index={0:"immature.B.cell--0",1:"immature.B.cell--1",2:"immature.B.cell--2"})
for cp in cp_list2:
    df=pd.DataFrame()
    df=pd.read_csv(f'/data1/qintian/science_repeat/data_mouse/GSE214390/result/sum_frame/{cp}.csv',header=0)
    df=df.rename(index={0:f'{cp}--0',1:f'{cp}--1',2:f'{cp}--2'})
    df0=pd.concat([df0,df],axis=0,join='inner')

df1 = pd.read_csv("/data1/qintian/science_repeat/data_mouse/GSE214390/result/nmr/immature.B.cell.csv",header=0)
df1 = df1.loc[[1]]
df1=df1.rename(index={1:"immature.B.cell--4"})
for cp in cp_list2:
    df=pd.DataFrame()
    df=pd.read_csv(f'/data1/qintian/science_repeat/data_mouse/GSE214390/result/nmr/{cp}.csv',header=0)
    df = df.loc[[1]]
    df=df.rename(index={1:f'{cp}--4'})
    df1=pd.concat([df1,df],axis=0,join='inner')

df2 = pd.read_csv("/data1/qintian/science_repeat/data_mouse/GSE214390/result/mouse/immature.B.cell.csv",header=0)
df2 = df2.loc[[1]]
df2=df2.rename(index={1:"immature.B.cell--4"})
for cp in cp_list2:
    df=pd.DataFrame()
    df=pd.read_csv(f'/data1/qintian/science_repeat/data_mouse/GSE214390/result/mouse/{cp}.csv',header=0)
    df = df.loc[[1]]
    df=df.rename(index={1:f'{cp}--4'})
    df2=pd.concat([df2,df],axis=0,join='inner')

gene_list=list(df0)
gene_list = [x for x in gene_list  if x != "aa"]

q,w,e,r=[],[],[],[]
for ge in gene_list:
    for cp in cp_list:
        if (exp[ge][cp]>0.01) and (df0[ge][f'{cp}--0']>0) and (df1[ge][f'{cp}--4']> 0) :
        #if (exp[ge][cp]>0.05) and (df0[ge][f'{cp}--0']>0) and (df1[ge][f'{cp}--4']> 0) :
            q.append(ge)
            w.append(cp)
            e.append(df0[ge][f'{cp}--0']/df1[ge][f'{cp}--4'])
            r.append(df0[ge][f'{cp}--1']/df1[ge][f'{cp}--4'])

n_df = pd.DataFrame({
    'gene': q,
    'celltype': w,
    'Dsp/V': e,
    'Dsam/V':r

    })

x_values = []
y_values = []

for i, row in n_df.iterrows():
    gene = row['gene']
    celltype = row['celltype']
    value = exp2[gene][celltype]


    x_values.append(value)
  
n_df['exp']=x_values

sorted_df1 = n_df.sort_values(by='Dsp/V', ascending=False)
sorted_df2 = n_df.sort_values(by='Dsam/V', ascending=False)
sorted_df1_len = len(sorted_df1)
sorted_df2_len = len(sorted_df2)

top_1_percent_rows = int(sorted_df1_len * 0.01)

sorted_df1_top_1_percent = sorted_df1.head(top_1_percent_rows)

bottom_5_percent_rows = int(len(sorted_df2) * 0.05)

sortdf0 = sorted_df1_top_1_percent[~sorted_df1_top_1_percent.isin(sorted_df2.head(bottom_5_percent_rows))].dropna()

sortdf0.to_csv(f'/data1/qintian/science_repeat/data_mouse/GSE214390/result/nmr/gogene_newmo.csv',index=False)
n_df.to_csv(f'/data1/qintian/science_repeat/data_mouse/GSE214390/result/nmr/gogene_newmo_allgenes.csv',index=False)

q,w,e,r=[],[],[],[]
for ge in gene_list:
    for cp in cp_list:
        if (exp2[ge][cp]>0.01) and (df0[ge][f'{cp}--0']>0) and (df2[ge][f'{cp}--4']> 0) :
        #if (exp[ge][cp]>0.05) and (df0[ge][f'{cp}--0']>0) and (df1[ge][f'{cp}--4']> 0) :
            q.append(ge)
            w.append(cp)
            e.append(df0[ge][f'{cp}--0']/df2[ge][f'{cp}--4'])
            r.append(df0[ge][f'{cp}--1']/df2[ge][f'{cp}--4'])

n_df2 = pd.DataFrame({
    'gene': q,
    'celltype': w,
    'Dsp/V': e,
    'Dsam/V':r

    })


sorted_df1 = n_df2.sort_values(by='Dsp/V', ascending=False)
sorted_df2 = n_df2.sort_values(by='Dsam/V', ascending=False)
sorted_df1_len = len(sorted_df1)
sorted_df2_len = len(sorted_df2)

top_1_percent_rows = int(sorted_df1_len * 0.01)

sorted_df1_top_1_percent = sorted_df1.head(top_1_percent_rows)

bottom_5_percent_rows = int(len(sorted_df2) * 0.05)

sortdf1 = sorted_df1_top_1_percent[~sorted_df1_top_1_percent.isin(sorted_df2.head(bottom_5_percent_rows))].dropna()

sortdf1.to_csv(f'/data1/qintian/science_repeat/data_mouse/GSE214390/result/mouse/gogene_newmo.csv',index=False)
n_df2.to_csv(f'/data1/qintian/science_repeat/data_mouse/GSE214390/result/mouse/gogene_newmo_allgenes.csv',index=False)