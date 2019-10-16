"""
Single Cell RNA Seq in Python (Jupyter Lab) 16th October 2019

Launch Jupyter Lab in folder that you want to work from (and in environment that you want to work in - might have to conda install
Juptyper Lab in each conda environment)

SCANPY: large scale single cell gene expression data analysis, Wolf et al (http://github.com/theislab/Scanpy)
https://icb-scanpy.readthedocs-hosted.com/en/stable/ -> API (preprocessing, tools and plotting are key functions)
Faster than Seurat and Cell Ranger R kit
Built on ANNDATA class - stores data matrix X, cell metadata info "observed", gene metadata info "variable"
Incorporate advanced machine learning techniques from sklearn site tensorflow
ANNDATA object, provided by AnnData class (annotated data):
	Data matrix X (raw data - counts or UMIs; numpy array, scipy sparse matrix) - observations/cells are rows, 
	features/variables are columns
	Gene metadata - dataframe-like object, is it mitochondrial or ercc?
	Cell metadata - dataframe-like object, cell IDs, batch information, cell cycle info, patient/sample info
	Unstructured info - dictionary-like
	Key difference with R - genes are rows, samples are columns. Have to transpose if going between R and Python

Case study used in tutorial - Lineage tracking reveals dynamic relationships of T cells in colorectal cancer (Zhange et al Nature 2018)
	11138 cells
	12 patients with CRC - 4 microsatellite instable (MSI), 8 microsatellite stable (MSS)
	20 T cell subsets
	Smart-Seq2 protocol
	Average 1.25 million unique read
	MSI/MSS response to checkpoint inhibition
"""
#lets you put plots in Jupyter Lab/Notebook along with code
%matplotlib inline

#anndata tutorial and colorectal cancer tutorial (iPython notebooks in Jupyter Lab)
#anndata tutorial and colorectal cancer tutorial as html files - type up code!
