import sys #this is for adding the module searching sites on my laptop, it can be changed or commented
sys.path.append('/Library/Frameworks/Python.framework/Versions/3.7/lib/python3.7/site-packages')#this is for adding the module searching sites on my laptop, it can be changed or commented
import mygene

#usage: python3 IDconvert.py <the input gene list>
#convert all id to gene symbol

#please input the genes list as a string as: geneID1\ngeneID2\n...geenIDX\n, which will be passed to this script as a argument

if len(sys.argv)>1 and len(sys.argv[1])>0:
	genelist=sys.argv[1].split('\\n')
	genelist=[i for i in genelist if i!='']
	mg=mygene.MyGeneInfo()
	result={}
	convert={}
	convert=mg.querymany(genelist,species='human',scopes='entrezgene, ensembl.gene, symbol, alias, refseq, unigene, accession, uniprot, pdb, interpro, pharmgkb, reporter',fields='symbol,name')
	if convert is not None:
		for dict in convert: #here the repeat genes are finally tiling to a uniq gene symbol
			if 'notfound' in dict:
				result[dict['query']]='not found' 
			else:
				result[dict['query']]=dict['symbol']+"\t"+dict['name']
	print(result)

else:
	print("The input contains no genes\n")