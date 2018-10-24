!/usr/bin/env python2
!/usr/bin/python
 -*- coding: utf-8 -*-
 Import libraries
import sys
import os
from collections import defaultdict
from collections import OrderedDict
import csv
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
"""
Created on Mon Oct  8 16:23:54 2018

@author: tomas

This script will plot histograms from distribution of uridylated sequences in all uridylated libraries.
"""
# Read in csv data and examine first 10 rows
TTT_distrib = pd.read_csv('Raw_%s.csv' % FILENAME, delimiter="\t")
TTT_distrib.head(10)
print("Plotting Histogram from results stored in file {}").format('Raw_%s.csv' % FILENAME)

## matplotlib histogram
# plt.hist(TTT_distrib['uridyl_length'], color = 'blue', edgecolor = 'black',
#          bins = int(180/5))
## Seaborn histogram
#sns.distplot(TTT_distrib['uridyl_length'], hist=True, kde=False, 
#             bins=int(100), color = 'blue',
#             hist_kws={'edgecolor':'black'})
## Shaded Density plot
#sns.distplot(subset['arr_delay'], hist = False, kde = True,
#                 kde_kws = {'shade': True, 'linewidth': 3}, 
#                  label = airline)
## Density Plot with Rug Plot
# WARNING! Rugplot will get stuck when too many items to plot and will use huge amount of memory
#sns.distplot(TTT_distrib['uridyl_length'], kde = True, rug = True,
#             color = 'darkblue', 
#             kde_kws={'linewidth': 3, "label": "FILENAME"},
#             rug_kws={'color': 'black'})
## Density Plot and Histogram
for line in open(file):
    sns.distplot(TTT_distrib['uridyl_length'], hist=True, kde=False, 
             bins=int(100), color = 'darkblue', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4})

## Add labels
plt.title('Histogram of Uridylation lengths', size = 18)
plt.xlabel('uridylation lengths', size = 16)
plt.ylabel('uridylation counts', size = 16)
