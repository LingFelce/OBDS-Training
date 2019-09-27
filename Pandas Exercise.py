#!/usr/bin/env python
# coding: utf-8

# # Pandas & Plotting Exercises

# Load read count dataset from week 1 and write it out as an Excel file
# 
# Moved featureCounts.summary from obds-training/exercises/rnaseq to obds-training/lingf/week1/rnaseq
# 
# Make sure move to right directory

# In[6]:


cd /ifs/obds-training/lingf/week1/rnaseq/


# In[74]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#makes plots display nicely in notebook
get_ipython().run_line_magic('matplotlib', 'inline')


# In[75]:


read_counts = pd.read_csv('featureCounts.summary', sep='\t', comment='#') #remove header at top, set column names 
read_counts.info
read_counts.to_csv('read_counts.csv')
read_counts.to_excel('read_counts.xlsx', sheet_name='Sheet1')


# Plot histogram and density plot of read counts per gene for each sample

# In[41]:


read_counts 
counts = read_counts.set_index('Geneid') #set Geneid as row name
counts = counts.iloc[:,-6:] #kept read counts for 6 samples (got rid of chromosome start/end, strand info etc)
counts


# In[61]:


d = counts["ERR1755082.bam"]
p = sns.distplot(d, color="m") #for histogram only kde=False, for density plot hist=False
plt.title("ERR1755082.bam", fontsize=20)
plt.xlabel("Counts per gene", fontsize=15)
plt.ylabel("Frequency",fontsize=15)


# In[118]:


for x in counts.columns: #to plot histograms for all samples
    if x.endswith('.bam'): #look for column heading which ends with .bam
        print(x)
        sns.distplot(counts[x], color="m")
        plt.xlabel("Counts per gene", fontsize=15)
        plt.ylabel("Frequency",fontsize=15)
        plt.show()        


# Plot scatter plot of correlation of replicates

# In[70]:


g = sns.pairplot(counts) #samples plotted against each other, look at scatter plot 
                         #2, 3 and 4 cluster tightly, 5, 6, and 7 cluster
                         #therefore 2, 3, and 4 are triplicate (CD4); 5, 6, and 7 are triplicate (CD8)


# Normalise read counts to total reads per sample

# In[83]:


def normalise_col(counts):
    return counts / counts.sum() #counts divided by sum of counts per row
norm_counts = counts.apply(normalise_col)
norm_counts


# Log transform normalised read counts

# In[86]:


log_counts = np.log2(norm_counts + 1) #add 1 to all values as some are 0, so can do log of 1
log_counts


# Create new dataframe of 10 genes with highest average expression for each condition

# In[100]:


CD4_counts = log_counts.iloc[:,0:3] #separate tables by condition
CD4_counts
CD8_counts = log_counts.iloc[:,3:7]
CD8_counts


# In[106]:


CD4_sorted_counts = CD4_counts.sort_values(by=['ERR1755082.bam'], ascending=False) #sort values by sample 1
top10_CD4_sorted_counts = CD4_sorted_counts.iloc[0:10,:] #select top 10
top10_CD4_sorted_counts


# In[107]:


CD8_sorted_counts = CD8_counts.sort_values(by=['ERR1755085.bam'], ascending=False) #sort values by sample 1
top10_CD8_sorted_counts = CD8_sorted_counts.iloc[0:10,:] #select top 10
top10_CD8_sorted_counts


# Identify top 100 most variable genes

# In[112]:


variance = log_counts.var(axis=1) #calculate variance across row for each gene
sorted_variance = variance.sort_values(ascending=False) #sort by descending
sorted_variance = sorted_variance.iloc[0:100] #select top 100
sorted_variance


# Plot line plot of expression of genes across all samples

# In[120]:


for x in log_counts.columns: #to plot histograms for all samples
    if x.endswith('.bam'): #look for column heading which ends with .bam
        print(x)
        sns.distplot(log_counts[x], hist=False, color='m')
        plt.xlabel("Log Counts", fontsize=15)
        plt.ylabel("Frequency",fontsize=15)
        plt.show() 


# Plot heatmap for top 100 most variable genes across all samples

# In[126]:


#merge gene list from sorted_variance with log_counts table
sorted_variance.name = 'variance' #name variance column in series
result = pd.merge(sorted_variance, log_counts, left_index=True, right_index=True) #merge series and table by left and right index
result = result.drop('variance', axis=1) #get rid of variance column to only show log counts and Geneid
result


# In[133]:


sns.heatmap(result, annot=False, cmap="Blues")


# Plot violin plot for top 5 most variable genes across two conditions

# In[168]:


top5 = result.iloc[0:5,:]
top5_CD4 = top5.iloc[:,0:3]
top5_CD8 = top5.iloc[:,3:6]
top5['CD4'] = top5_CD4.mean(axis=1)
top5['CD8'] = top5_CD8.mean(axis=1)
top5 = top5.iloc[:,6:8]
top5


# In[173]:


sns.violinplot(data=top5, palette="Set3", bw=.2, cut=1, linewidth=1) #combined violin plot


# In[ ]:




