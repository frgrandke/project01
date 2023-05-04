import json

# uses the modified graviton.py script in this folder!
from graviton import *
import numpy
import pandas as pd
import zipfile

# Obtain a session

degs  = pd.read_csv("results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep= "\t")
degs = degs[(degs["contrast"] == "ADvsHC") | (degs["contrast"] == "PDvsHC") | (degs["contrast"] == "PDMCIvsHC") | (degs["contrast"] == "MCIvsHC")]
print(degs)

def combinePathways(filename):
	dfs = pd.DataFrame()
	files = ["GO_-_Biological_Process.txt", "GO_-_Cellular_Component.txt", "GO_-_Molecular_Function.txt", "KEGG_-_Pathways.txt", "Pfam_-_Protein_families.txt", "Reactome_-_Pathways.txt"]
	zf = zipfile.ZipFile(filename) 

	for curr_file in files: 
		try: 
			df = pd.read_csv(zf.open(curr_file), sep="\t")
		except: 
			continue

		df = df[df["P-value"] < 0.05]
		df["Category"] = curr_file

		if (dfs.size == 0): 
			dfs = df
		else: 
			dfs = pd.concat([dfs, df])
	return dfs

def getPathwayAnalysisGSEA(degs):
	genes_up["-log10_p_adj"] = -numpy.log10(genes_up["p_adj.loc"])
	genes_up["abs_logFC"] = abs(genes_up["logFC"])

	genes_up.sort_values(by=['-log10_p_adj', 'abs_logFC'], ascending=False)

	numpy.savetxt(r'genes_up.txt', genes_up["gene"], fmt='%s', delimiter = "\t")

	if genes_up["gene"].size <=2:
		return pd.DataFrame()
	try : 
		scores = uploadFile(key, "genes_up.txt", organism)["id"]
		setupEnrichment(key, 'gsea', scores, mrnaCategories)
		result = (runJob(key)['enrichment']['id'])
	except: 
		return None

	downloadResult(key, result, 'mrnaAllSamples.gsea.zip')

	dfs = combinePathways("mrnaAllSamples.gsea.zip")
	
	return dfs



type = "ORA"
organism = "9606"

key = ""

if key == "":
    print("Getting session key...")
    key = getSession(organism)
    print("Key: {}".format(key))

mrnaCategories = ['9606-gene-go-biologicalprocess',
  '9606-gene-go-cellularcomponent',
  '9606-gene-go-molecularfunction', 
  '9606-gene-kegg-pathways',
  '9606-gene-reactome-pathways',
  '9606-gene-pfam-proteinfamilies', 
  '9606-gene-wikipathways']

print(list(degs.columns))

celltypes = list(set(degs["cluster_id"]))
comparisons = list(set(degs["contrast"]))

results_path = "results/genetrail3/"

for celltype in celltypes: 
	for comparison in comparisons: 
		print(" ")
		print(celltype)
		print(comparison)
		print(" ")
		genes_up = degs[(degs["cluster_id"] == celltype) & (degs["contrast"] == comparison)] 
		dfs = getPathwayAnalysisGSEA(genes_up)
		
		directory = "/".join([comparison])
		path = os.path.join(results_path, directory)
		if not os.path.isdir(path): 
			os.mkdir(path)
		path = os.path.join(path, celltype)
		if not os.path.isdir(path): 
			os.mkdir(path)
		if dfs is not None: 
			dfs.to_csv(os.path.join(path, "patways_upregulated_genes.csv"))
		
		genes_up = degs[(degs["cluster_id"] == celltype) & (degs["contrast"] == comparison)& (degs["logFC"] < 0)] 
		dfs = getPathwayAnalysisGSEA(genes_up)
		directory = "/".join([comparison])
		path = os.path.join(results_path, directory)
		if not os.path.isdir(path): 
			os.mkdir(path)
		path = os.path.join(path, celltype)
		if not os.path.isdir(path): 
			os.mkdir(path)
		if dfs is not None: 
			dfs.to_csv(os.path.join(path, "patways_downregulated_genes.csv"))
		
		genes_up = degs[(degs["cluster_id"] == celltype) & (degs["contrast"] == comparison) & (degs["logFC"] > 0)]
		dfs = getPathwayAnalysisGSEA(genes_up)
		directory = "/".join([comparison])
		path = os.path.join(results_path, directory)
		if not os.path.isdir(path): 
			os.mkdir(path)
		path = os.path.join(path, celltype)
		if not os.path.isdir(path): 
			os.mkdir(path)
		if dfs is not None: 
			dfs.to_csv(os.path.join(path, "patways_deregulated_genes.csv"))


