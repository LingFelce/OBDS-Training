#11_Data science in Python 27th September 2019

"""
Data science - analysis and visualisation of data
Pandas is Python package for data analysis, built using numerical analysis package NumPy
Plotting in Python - many packages built on Matplotlib; Seaborn is much nicer - http://seaborn.pydata.org/examples/index.html
Use Jupyter Notebook (previously iPython Notebook) - launch from Terminal, more similar to creating .Rmd files,
can integrate output with input and export as HTML with data plots etc. Good for producing reports
Jupyter Notebook will be replaced by JupyterLab
Pandas - converts data structures into DataFrame objects; insert/delete columns, indexing, slicing and subsetting

"""

#NumPy - array, numerical datatypes (int64, float64),mathematical functions, slicing and indexing
import numpy as np
x = np.array([[1,2],[3,4]], dtype = np.float64)
b = x[:2,1:3] #slicing

#Pandas series - one-dimensional vector
import pandas as pd
s = pd.Series([1,2,3,4,5,6,7,8]) #pandas creates default integer index
s = pd.Series(np.random.randn(5),
	index=['a','b','c','d','e']) #add own index
s.index #attributes of Series class
s.name

#Pandas dataframe - two-dimensional, includes axis labels, any Python object, different data type per column
#dictionary of series/dictionary of dictionary, list of dictionary/dictionary of lists, dictionary of tuples
d = {'one': pd.Series([1., 2., 3.], #Dataframe from dictionary of series
	index=['a', 'b', 'c']),
	'two': pd.Series([1., 2., 3., 4.], 
	index=['a', 'b', 'c', 'd'])}
df= pd.DataFrame(d)

# attributes
df.index #row labels
df.columns #column labels
df.shape #dataframe size

df.head() #Viewing dataframes
df.tail()

df.info() #Summary

# Adding columns
df['three'] = df['one'] * df['two']
df['flag'] = df['one'] > 2

# Assign –like dplyrmutate
iris.assign(sepal_ratio= 
	iris['SepalWidth'] / 
	iris['SepalLength'])

# Deleting columns
del df['two']
three = df.pop('three')

df[col] #select column - series
df.loc[label] #select row by label - series
df.iloc[loc] #select row by integer location - series
df[5:10] #slice rows - DataFrame
df[bool_vec] #select rows by boolean vector - DataFrame

df[‘col1’] #Select single column (don't have to do print)
df[[‘col2’, ‘col1’]] #Select multiple columns (list)

df[:3]df.loc[‘rowA’,:] #Slice data frame rows / cols
df.loc[‘row3':, ‘col2':‘col4']
df.iloc[2:]
df1.iloc[1:5, 2:4]

#Boolean indexing
df[df.A> 0] #select rows on column values
df[df> 0] #select values on condition

df["biotype"] = df["raw_biotype"].astype("category") #categorical variables, can be ordered 
df["biotype"].cat.categories= ["protein_coding", "rRNA", "lincRNA"]


df.dropna(how='any') #drop rows with nan
df.fillna(value=5) #Replace nan with specific value
pd.isna(df) #Create Boolean mask

df* 5 #Element-wise mathematical operations
df**3

df1 & df2 #Boolean logic

np.log(df) #Numpyfunctions

df.describe() #Statistical summary

#statistical functions
df.mean() #per column

df.mean(1) #per row

# Matrix algebra
df.T #transpose of matrix
df.dot(df) #matrix multiplication

#apply functions along axis of data frame
df.apply(np.sqrt) #apply to entire dataframe
df.apply(np.sum, axis=1) #apply to cols
df.apply(np.sum, axis=0) #apply to rows

df.sort_values(by=['col1']) #Sort by single column
df.sort_values(by=['col1', 'col2']) #Sort by multiple columns
df.sort_values(by='col1', ascending=False) #Descending order
df.sort_values(by='col1', ascending=False, na_position='first') #NAs first

#combining datasets - Linux like concatenation - concat
#SQL style merges - join() and merge() - join is by position, merge is by rownames/colnames?
#merge left - keys from left frame only; right - keys from right frame only; outer - union of keys from both frames
#inner - intersection of keys from both frames

pd.concat([df[:3], df[3:7], df[7:]]) #Concatenate rows
result = pd.merge(df1, df2, on='key') #Merge (inner by default)
result = pd.merge(df1, df2, on='key', 
		how='left')
result = df1.join(df2, how='outer') #Join on index
result = df1.join(df2, on='key') #Join on key

#Group rows by categorical variable in column A then calculate the sum for each group
#Produces a dataframe where rows are groups in A and columns are sum of each numerical column
df.groupby('A').sum()

melted = df.melt(id_vars=['country'], var_name='year', value_name='cases') #tidy data
df= melted.pivot(index='country', columns='year') #Unmelt

pd.read_csv('foo.csv') #Read and write to CSV
df.to_csv('foo.csv')
pd.read_csv('foo.tsv', sep='\t‘, header=0) #tsv
pd.read_table(fpath) #tab default
pd.read_excel('foo.xlsx', 'Sheet1', index_col=None, na_values=['NA']) #Excel
df.to_excel('foo.xlsx',sheet_name='Sheet1')
pd.read_hdf('foo.h5', 'df') #HDF5
df.to_hdf('foo.h5', 'df')

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipyimport stats

x = np.random.normal(size=100)
sns.distplot(x)

sns.distplot(x, bins=20, kde=False, rug=True)  #Select number of bins, remove kernel density estimation, add rug plot

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="darkgrid")
tips = sns.load_dataset("tips")
sns.relplot(x="total_bill", y="tip", hue="smoker", data=tips) #scatter plot example
sns.relplot(x="total_bill", y="tip", hue="smoker", style="time", data=tips)
sns.relplot(x="total_bill", y="tip", hue="smoker", col="time", data=tips) #multiple relationship with Facets

iris = sns.load_dataset("iris") #pairs plot
g = sns.PairGrid(iris, hue="species")
g.map_diag(plt.hist)
g.map_offdiag(plt.scatter)
g.add_legend();

fmri= sns.load_dataset("fmri")
sns.relplot(x="timepoint", y="signal", hue="event", kind="line", data=fmri) #lineplot
sns.relplot(x="timepoint", y="signal", hue="region", style="event", dashes=False, markers=True, kind="line", data=fmri)
sns.relplot(x="timepoint", y="signal", hue="subject", col="region", row="event", height=3, kind="line", estimator=None, data=fmri)

g = sns.catplot(x="day", y="total_bill", kind="violin", inner=None, data=tips) #combined violin and swarm plot
sns.swarmplot(x="day", y="total_bill", color="k", size=3, data=tips, ax=g.ax)

#HEATMAP/CLUSTERMAP
import pandas as pd
import seaborn as sns 
sns.set()
# Load the brain networks dataset
df= sns.load_dataset("brain_networks", header=[0, 1, 2], index_col=0)
# Select a subset of the networks
used_networks= [1, 5, 6, 7, 8, 12, 13, 17]
used_columns= (df.columns.get_level_values("network").astype(int).isin(used_networks))
df= df.loc[:, used_columns]
# Create a categorical palette to identify the networks
network_pal= sns.husl_palette(8, s=.45)
network_lut= dict(zip(map(str, used_networks), network_pal))
# Convert the palette to vectors that will be drawn on the side of the matrix
networks = df.columns.get_level_values("network")
network_colors= pd.Series(networks, index=df.columns).map(network_lut)
# Draw the full plot
sns.clustermap(df.corr(), center=0, cmap="vlag", row_colors=network_colors, col_colors=network_colors, linewidths=.75, figsize=(13, 13))

