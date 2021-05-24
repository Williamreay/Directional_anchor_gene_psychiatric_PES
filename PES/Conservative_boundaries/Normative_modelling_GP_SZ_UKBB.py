#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  3 11:25:18 2021

@author: williamreay

"""

## Normative modelling script

import pandas as pd
import numpy as np
from pynm import pynm
import seaborn as sns
import matplotlib.pyplot as plt

import patsy as pat
import statsmodels.api as sm
import time

np.random.seed(666)

## FES network PES ##

FES_df=pd.read_csv("~/Desktop/SZ_PES_mtCOJO_norm/PES/UKBB_SZ/Norm_modelling/FES_norm_modelling_input.csv")

## Initialise with data - define outcome and coefficients

py_nm = pynm.PyNM(FES_df,'score','group',
            conf = 'PRS_0_5',                           #PRS main variable for LOESS/centile
            confounds = ['Age','C(Sex)', 'PC1','PC2', 'PC3', 'PC4', 'PC5']) #multivarite confounds for GP model


#I# nitialize bins for LOESS and Centiles models
x, xcount, z, zm, zstd, zci = py_nm.create_bins()

## Run LOESS model
loess = py_nm.loess_normative_model()

## Run centiles model

centiles = py_nm.centiles_normative_model()

## Run Gaussian process

gp = py_nm.gp_normative_model() 

## Derive output

data = py_nm.data


## Plot GP only with PRS as term


data.sort_values('PRS_0_5',inplace=True)

plt.figure(figsize=(15,8))
sns.scatterplot(data=data,x='PRS_0_5', y='score', label='Data',hue='group')
plt.plot(data.PRS_0_5, data.GP_nmodel_pred, color='black', label='Prediction',lw=2)
plt.fill_between(np.squeeze(data.PRS_0_5),np.squeeze(data.GP_nmodel_pred) - 2*data.GP_nmodel_sigma,np.squeeze(data.GP_nmodel_pred) + 2*data.GP_nmodel_sigma,alpha=.2, fc='grey', ec='None', label='95% CI')

plt.xlabel('$PRS$')
plt.ylabel('$PES$')
plt.legend(loc='upper right')
plt.title('FES network PES vs SZ PRS')
plt.xlim(-0.006886405,-0.006738266)




# plt.figure(figsize=(15,8))

# ##Plot centiles
# plt.fill_between(x,z[:, 5],z[:, 95],facecolor='grey',alpha=0.2,label='5-95')#,zorder=1)
# plt.fill_between(x,z[:, 25],z[:, 75],facecolor='grey',alpha=0.4,label='25-75')#,zorder=2)

# #Plot participants
# sns.scatterplot(x = 'PRS_0_5',y='score',data=data, hue='group',style='group')

# plt.plot(x, z[:, 50], '-', color='black',label='nmodel')

# #Plot bins with less than 20 participants
# for bin_sel in np.argwhere((xcount <= 100) + (xcount <= 100)).flatten():
#     plt.plot(np.array([x[bin_sel], x[bin_sel]]),
#                 [z.min() - 1, z.max() + 1], '--k', alpha=0.5,color='red')

# plt.ylabel('PES')
# plt.xlabel('PRS')
# plt.xlim(-0.006886405,-0.006738266)
# plt.legend(loc='lower left')
# plt.show()


## GP corrected for other confounders

# plt.figure(figsize=(15,8))
# sns.scatterplot(data=data,x='PRS_0_5', y='GP_nmodel_residuals',hue='group',style='group')

# plt.xlabel('$PRS$')
# plt.ylabel('$PES$')
# plt.legend(loc='upper right')
# plt.title('Corrected')

# plt.figure(figsize=(15,8))
# sns.scatterplot(data=data,x='PRS_0_5', y='score',hue='group',style='group')

# plt.xlabel('$PRS$')
# plt.ylabel('$PES$')
# plt.legend(loc='upper right')
# plt.title('Uncorrected')


