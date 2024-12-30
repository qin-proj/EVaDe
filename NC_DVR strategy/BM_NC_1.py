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
#Rpm = r(Vsp ~ Vsam), Rps = r(Vsp ~ Vres), Rms = (Vsam ~ Vres)
#NMR
ge_li,Rpm,Rps,Rms,pv1,pv2,pv3,expr=[],[],[],[],[],[],[],[]
for ge in gene_list:
    x,y,z=[],[],[]
    for cp in cp_list:
        #if (df0[ge][f'{cp}--0'] > 0) and (df1[ge][f'{cp}--4']>0):
        if (df0[ge][f'{cp}--0'] > 0) and (df1[ge][f'{cp}--4']>0) and (exp[ge][cp]>0):
            x.append(df0[ge][f'{cp}--0']/exp[ge][cp])
            y.append(df0[ge][f'{cp}--1']/exp[ge][cp])
            z.append(df1[ge][f'{cp}--4']/exp[ge][cp])
    if len(x)>5:
        slope, intercept, r1, p1, std_err = stats.linregress(x, y)
        slope, intercept, r2, p2, std_err = stats.linregress(x, z)
        slope, intercept, r3, p3, std_err = stats.linregress(y, z)
        ge_li.append(ge)
        Rpm.append(r1)
        pv1.append(p1)
        Rps.append(r2)
        pv2.append(p2)
        Rms.append(r3)
        pv3.append(p3)

corr_df = pd.DataFrame({
    'gene': ge_li,
    'Rpm': Rpm,
    'p1':pv1,
    'Rps': Rps,
    'p2':pv2,
    'Rms': Rms,
    'p3':pv3,
    })
corr_df.to_csv(f'/data1/qintian/science_repeat/data_mouse/GSE214390/result/nmr/gogene_all.csv',index=False)

#mouse
ge_li,Rpm,Rps,Rms,pv1,pv2,pv3,expr=[],[],[],[],[],[],[],[]
for ge in gene_list:
    x,y,z=[],[],[]
    for cp in cp_list:
        #if (df0[ge][f'{cp}--0'] > 0) and (df1[ge][f'{cp}--4']>0):
        if (df0[ge][f'{cp}--0'] > 0) and (df2[ge][f'{cp}--4']>0) and (exp2[ge][cp]>0):
            x.append(df0[ge][f'{cp}--0']/exp2[ge][cp])
            y.append(df0[ge][f'{cp}--1']/exp2[ge][cp])
            z.append(df2[ge][f'{cp}--4']/exp2[ge][cp])
    if len(x)>5:
        slope, intercept, r1, p1, std_err = stats.linregress(x, y)
        slope, intercept, r2, p2, std_err = stats.linregress(x, z)
        slope, intercept, r3, p3, std_err = stats.linregress(y, z)
        ge_li.append(ge)
        Rpm.append(r1)
        pv1.append(p1)
        Rps.append(r2)
        pv2.append(p2)
        Rms.append(r3)
        pv3.append(p3)

corr_df = pd.DataFrame({
    'gene': ge_li,
    'Rpm': Rpm,
    'p1':pv1,
    'Rps': Rps,
    'p2':pv2,
    'Rms': Rms,
    'p3':pv3,
    })
corr_df.to_csv(f'/data1/qintian/science_repeat/data_mouse/GSE214390/result/mouse/gogene_all.csv',index=False)