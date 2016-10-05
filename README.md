# CREEDS: CRowd Extracted Expression of Differential Signatures
The web portal for serving the crowdsourced gene expression signatures from Gene Expression Omnibus: http://amp.pharm.mssm.edu/creeds/. The web portal allows users to query, download, and visualize these gene expression signatures.

## Abstract
Gene expression data are accumulating exponentially in public repositories. Reanalysis and integration of themed collections from these studies may provide new insights, but requires further human curation. Here we report a crowdsourcing project to annotate and reanalyse a large number of gene expression profiles from [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/). Through a massive open online course on Coursera, over 70 participants from over 25 countries identify and annotate 2,460 single-gene perturbation signatures, 839 disease versus normal signatures, and 906 drug perturbation signatures. All these signatures are unique and are manually validated for quality. Global analysis of these signatures confirms known associations and identifies novel associations between genes, diseases and drugs. The manually curated signatures are used as a training set to develop classifiers for extracting similar signatures from the entire GEO repository. We develop a web portal to serve these signatures for query, download and visualization.

## Publication
Wang, Z., Monteiro, C. D., Jagodnik, K. M., Fernandez, N. F., Gundersen, G. W., ... & Ma'ayan, A. (2016) **Extraction and Analysis of Signatures from the Gene Expression Omnibus by the Crowd.**    
_Nature Communications_ doi: [10.1038/ncomms12846](http://dx.doi.org/10.1038/ncomms12846) PMID: [27667448](https://www.ncbi.nlm.nih.gov/pubmed/27667448)

## API Documentation
Three API endpoints have been implemented searching, querying, and retrieving associated metadata for gene expression signatures. For more details, visit http://amp.pharm.mssm.edu/CREEDS/#help/#api-doc.
Tutorials walking through all the API endpoints with examples are also available in [Python](http://nbviewer.jupyter.org/github/maayanlab/creeds/blob/master/Example_API_usage.ipynb) and [R](http://rpubs.com/wangz10/177826).

The site can also be accessed through http://amp.pharm.mssm.edu/creeds2/.

## In the News
- [CROWDSOURCING FOR SCIENTIFIC DISCOVERY: MOUNT SINAI RESEARCHERS FIND NOVEL WAYS TO ANALYZE DATA FOR DRUG AND TARGET DISCOVERY](http://www.newswise.com/articles/crowdsourcing-for-scientific-discovery-mount-sinai-researchers-find-novel-ways-to-analyze-data-for-drug-and-target-discovery)
- [Mount Sinai Crowdsources Analysis and Annotation of Gene Eexpression Profiles](http://www.clinicalomics.com/articles/mount-sinai-crowdsources-analysis-and-annotation-of-gene-expression-profiles/766)
