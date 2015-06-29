## R script to perform HC on data and produce a edgelist of the dendrogram

library(ape)


hc2network <- function(data, dist.method, linkage.method) {
	## data: a numeric matrix or data frame
	## dist.method should be one of ["euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"]
	## linkage.method should be one of ["ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)]
	if (dist.method == 'cosine') {
		library(MKmisc)
		D <- MKmisc::corDist(data, method = 'cosine')
	} else {
		D <- dist(data, method = dist.method)	
	}
	hc <- hclust(D, method = linkage.method)
	phylo.tree <- ape::as.phylo(hc)
	phylo.tree$edge # matrix
}
