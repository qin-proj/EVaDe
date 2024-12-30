import pyreadr
import pandas as pd
import scipy
from scipy import stats
import numpy as np
from matplotlib import pyplot
cp_list=["L2-3 IT", "L3-5 IT-1", "L3-5 IT-2", "L3-5 IT-3",  "L5-6 NP", "L6 CT", "L6 IT-1", "L6 IT-2", "L6B","LAMP5 LHX6", "LAMP5 RELN", "VIP", "ADARB2 KCNG1", "SST", "PVALB", "PVALB ChC","Astro","OPC","Oligo","Micro","Immune","Endo","PC","SMC","VLMC"]
cp_list2=["L3-5 IT-1", "L3-5 IT-2", "L3-5 IT-3",  "L5-6 NP", "L6 CT", "L6 IT-1", "L6 IT-2", "L6B","LAMP5 LHX6", "LAMP5 RELN", "VIP", "ADARB2 KCNG1", "SST", "PVALB", "PVALB ChC","Astro","OPC","Oligo","Micro","Immune","Endo","PC","SMC","VLMC"]

exp=pd.read_csv(f'/data1/qintian/science_repeat/ProcessedData/result0/anova/sum_frame_2sp_2/Human2/expression.csv',index_col=0)
exp2=pd.read_csv(f'/data1/qintian/science_repeat/ProcessedData/result0/anova/sum_frame_2sp_2/expression_Rhesus.csv',index_col=0)

df0 = pd.read_csv("/data1/qintian/science_repeat/ProcessedData/result0/anova/sum_frame_rh/L2-3 IT.csv",header=0)
#0:sp 1:sam 2:Residuals
df0=df0.rename(index={0:"L2-3 IT--0",1:"L2-3 IT--1",2:"L2-3 IT--2"})
for cp in cp_list2:
    df=pd.DataFrame()
    df=pd.read_csv(f'/data1/qintian/science_repeat/ProcessedData/result0/anova/sum_frame_rh/{cp}.csv',header=0)
    df=df.rename(index={0:f'{cp}--0',1:f'{cp}--1',2:f'{cp}--2'})
    df0=pd.concat([df0,df],axis=0,join='inner')

df1 = pd.read_csv("/data1/qintian/science_repeat/ProcessedData/result0/anova/sum_frame_rh/Human/L2-3 IT.csv",header=0)
df1 = df1.loc[[1]]
df1=df1.rename(index={1:"L2-3 IT--4"})
for cp in cp_list2:
    df=pd.DataFrame()
    df=pd.read_csv(f'/data1/qintian/science_repeat/ProcessedData/result0/anova/sum_frame_rh/Human/{cp}.csv',header=0)
    df = df.loc[[1]]
    df=df.rename(index={1:f'{cp}--4'})
    df1=pd.concat([df1,df],axis=0,join='inner')

df2 = pd.read_csv("/data1/qintian/science_repeat/ProcessedData/result0/anova/sum_frame_rh/Rhesus/L2-3 IT.csv",header=0)
df2 = df2.loc[[1]]
df2=df2.rename(index={1:"L2-3 IT--4"})
for cp in cp_list2:
    df=pd.DataFrame()
    df=pd.read_csv(f'/data1/qintian/science_repeat/ProcessedData/result0/anova/sum_frame_rh/Rhesus/{cp}.csv',header=0)
    df = df.loc[[1]]
    df=df.rename(index={1:f'{cp}--4'})
    df2=pd.concat([df2,df],axis=0,join='inner')

gene_list=list(df0)
gene_list = [x for x in gene_list  if x != "aa"]
#Human
q,w,e,r=[],[],[],[]
for ge in gene_list:
    for cp in cp_list:
        if (exp[ge][cp]>0.01) and (df0[ge][f'{cp}--0']>0) and (df1[ge][f'{cp}--4']> 0) :
        #if (exp[ge][cp]>0.2) and (df0[ge][f'{cp}--0']>0) and (df1[ge][f'{cp}--4']> 0) :
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

#add expression data
for i, row in n_df.iterrows():
    gene = row['gene']
    celltype = row['celltype']
    value = exp[gene][celltype]
    value2 = exp2[gene][celltype]

    x_values.append(value)
    y_values.append(value2)
n_df['exp']=x_values
n_df['exp2']=y_values

n_df.to_csv(f'/data1/qintian/science_repeat/ProcessedData/result0/anova/sum_frame_rh/Human/gogene_ratio_all.csv',index=False)
sorted_df1 = n_df.sort_values(by='Dsp/V', ascending=False)
sorted_df2 = n_df.sort_values(by='Dsam/V', ascending=False)
sorted_df1_len = len(sorted_df1)
sorted_df2_len = len(sorted_df2)

top_1_percent_rows = int(sorted_df1_len * 0.005)

sorted_df1_top_1_percent = sorted_df1.head(top_1_percent_rows)

bottom_5_percent_rows = int(len(sorted_df2) * 0.05)

sortdf0 = sorted_df1_top_1_percent[~sorted_df1_top_1_percent.isin(sorted_df2.head(bottom_5_percent_rows))].dropna()
sortdf0.to_csv(f'/data1/qintian/science_repeat/ProcessedData/result0/anova/sum_frame_rh/Human/gogene_newmo.csv.csv',index=False)
#Macaque
q,w,e,r=[],[],[],[]
for ge in gene_list:
    for cp in cp_list:
        if (exp2[ge][cp]>0.01) and (df0[ge][f'{cp}--0']>0) and (df2[ge][f'{cp}--4']> 0) :
        #if (exp2[ge][cp]>0.2) and (df0[ge][f'{cp}--0']>0) and (df2[ge][f'{cp}--4']> 0) :
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
n_df2.to_csv(f'/data1/qintian/science_repeat/ProcessedData/result0/anova/sum_frame_rh/Rhesus/gogene_ratio_all.csv',index=False)
sorted_df1 = n_df2.sort_values(by='Dsp/V', ascending=False)
sorted_df2 = n_df2.sort_values(by='Dsam/V', ascending=False)
sorted_df1_len = len(sorted_df1)
sorted_df2_len = len(sorted_df2)

top_1_percent_rows = int(sorted_df1_len * 0.005)

sorted_df1_top_1_percent = sorted_df1.head(top_1_percent_rows)

bottom_5_percent_rows = int(len(sorted_df2) * 0.05)

sortdf1 = sorted_df1_top_1_percent[~sorted_df1_top_1_percent.isin(sorted_df2.head(bottom_5_percent_rows))].dropna()
sortdf1.to_csv(f'/data1/qintian/science_repeat/ProcessedData/result0/anova/sum_frame_rh/Rhesus/gogene_newmo.csv.csv',index=False)
