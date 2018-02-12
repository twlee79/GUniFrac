
GUniFrac <- function (otu.tab, tree, alpha = c(0, 0.5, 1)) {
	# Calculate Generalized UniFrac distances. Unweighted and 
	# Variance-adjusted UniFrac distances will also be returned.
	#	
	# Args:
	#		otu.tab: OTU count table, row - n sample, column - q OTU
	#		tree: rooted phylogenetic tree of R class "phylo"
	#		alpha: parameter controlling weight on abundant lineages
	#
	# Returns:
	# 	unifracs: three dimensional array containing the generalized 
	#							UniFrac distances, unweighted UniFrac distance and 
	#							variance adjusted UniFrac distances. 
	#

	if (!is.rooted(tree)) stop("Rooted phylogenetic tree required!")
	
	# Convert into proportions
	otu.tab <- as.matrix(otu.tab)
	row.sum <- rowSums(otu.tab)
	otu.tab <- otu.tab / row.sum
	n <- nrow(otu.tab)
	
	# Construct the returning array
	if (is.null(rownames(otu.tab))) {
		rownames(otu.tab) <- paste("comm", 1:n, sep="_")
	}
	# d_UW: unweighted UniFrac, d_VAW: weighted UniFrac
	dimname3 <- c(paste("d", alpha, sep="_"), "d_UW", "d_VAW")
	unifracs <- array(NA, c(n, n, length(alpha) + 2),
				  dimnames=list(rownames(otu.tab), rownames(otu.tab), dimname3))
	for (i in 1:(length(alpha)+2)){
		for (j in 1:n){
			unifracs[j, j, i] <- 0
		}
	}	
	
	# Check OTU name consistency
	if (sum(!(colnames(otu.tab) %in% tree$tip.label)) != 0) {
		stop("The OTU table contains unknown OTUs! OTU names
					in the OTU table and the tree should match!" )
	}
	
	# Get the subtree if tree contains more OTUs
	absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
	if (length(absent) != 0) {
		tree <- drop.tip(tree, absent)
		warning("The tree has more OTU than the OTU table!")
	}
	
	# Reorder the otu.tab matrix if the OTU orders are different
	tip.label <- tree$tip.label
	otu.tab <- otu.tab[, tip.label]
	
	ntip <- length(tip.label)
	nbr <- nrow(tree$edge)	
	edge <- tree$edge
	edge2 <- edge[, 2]
	br.len <- tree$edge.length
	
	#  Accumulate OTU proportions up the tree	
	cum <- matrix(0, nbr, n)							# Branch abundance matrix
	for (i in 1:ntip) {
		tip.loc <- which(edge2 == i)
		cum[tip.loc, ] <- cum[tip.loc, ] + otu.tab[, i]	
		node <- edge[tip.loc, 1]						# Assume the direction of edge 
		node.loc <- which(edge2 == node)
		while (length(node.loc)) {
			cum[node.loc, ] <- cum[node.loc, ] + otu.tab[, i]		
			node <- edge[node.loc, 1]
			node.loc <- which(edge2 == node)
		}
	}
	
	# Calculate various UniFrac distances
	cum.ct <- round(t(t(cum) * row.sum)) 	# For VAW
	for (i in 2:n) {
		for (j in 1:(i-1)) {
			cum1 <- cum[, i]
			cum2 <- cum[, j]
			ind <- (cum1 + cum2) != 0
			cum1 <- cum1[ind]
			cum2 <- cum2[ind]		
			br.len2 <- br.len[ind]			
			mi <- cum.ct[ind, i] + cum.ct[ind, j]
			mt <- row.sum[i] + row.sum[j]			
			diff <- abs(cum1 - cum2) / (cum1 + cum2)		
			
			# Generalized UniFrac distance
			for(k in 1:length(alpha)){
				w <- br.len2 * (cum1 + cum2)^alpha[k]
				unifracs[i, j, k] <- unifracs[j, i, k] <- sum(diff * w) / sum(w)
			}			
			
			#	Variance Adjusted UniFrac Distance
			ind2 <- (mt != mi)
			w <- br.len2 * (cum1 + cum2) / sqrt(mi * (mt - mi))
			unifracs[i, j, (k + 2)] <- unifracs[j, i, (k + 2)] <- 
					sum(diff[ind2] * w[ind2]) / sum(w[ind2])		
			
			#	Unweighted UniFrac Distance
			cum1 <- (cum1 != 0)
			cum2 <- (cum2 != 0)			
			unifracs[i, j, (k + 1)] <- unifracs[j, i, (k + 1)] <- 
					sum(abs(cum1 - cum2) / (cum1 + cum2) * br.len2) / sum(br.len2)
		}
	}
	return(list(unifracs=unifracs))
}

 
PermanovaG <- function(formula, dat=NULL, ...) {
	# Distance based statistical test by combining multiple distance
	# matrices based on PERMANOVA procedure by taking the maximum of
	# pseudo-F statistics
	#	
	# Args:
	#		formula: left side of the formula is a three dimensional array
	#				 of the supplied distance matrices as produced by GUniFrac,
	#                Or a list of distance matrices.
	#		dat: data.frame containing the covariates
	#		...: Parameter passing to "adonis" function
	#
	# Returns:
	# 	p.tab (data.frame): rows - Covariates (covariate of interest should appear at the last)
	#        columns - F.model, p.values for individual distance matrices and the omnibus test, 
	#   aov.tab.list:  a list of the aov.tab from individual matrices 
	
	save.seed <- get(".Random.seed", .GlobalEnv)
	lhs <- formula[[2]]
	lhs <- eval(lhs, dat, parent.frame())
	rhs <- as.character(formula)[3]
	
	array.flag <- is.array(lhs)
	list.flag <- is.list(lhs)
	
	if (!array.flag & !list.flag) {
		stop('The left side of the formula should be either a list or an array of distance matrices!\n')
	}
	if (array.flag) {
		if (length(dim(lhs)) != 3) {
			stop('The array should be three-dimensional!\n')
		} else {
			len <- dim(lhs)[3]
		}
	}
	
	if (list.flag) {
		len <- length(lhs)
	}
	
	p.perms <- list()
	p.obs <- list()
	aov.tab.list <- list()
	for (i_ in 1:len) {
		assign(".Random.seed", save.seed, .GlobalEnv)
		if (array.flag) {
			Y <- as.dist(lhs[, , i_])
		}
		if (list.flag) {
			Y <- as.dist(lhs[[i_]])
		}

		obj <- adonis(as.formula(paste("Y", "~", rhs)), dat, ...)
		perm.mat <- obj$f.perms
		p.perms[[i_]] <- 1 - (apply(perm.mat, 2, rank) - 1) / nrow(perm.mat)
		p.obs[[i_]] <- obj$aov.tab[1:ncol(perm.mat), "Pr(>F)"]
		aov.tab.list[[i_]] <- obj$aov.tab
	}
	
	omni.pv <- NULL
	indiv.pv <- NULL
	for (j_ in 1:ncol(perm.mat)) {
		p.perms.j <- sapply(p.perms, function (x) x[, j_])
		p.obj.j <- sapply(p.obs, function (x) x[j_])
		omni.pv <- c(omni.pv, mean(c(rowMins(p.perms.j ) <= min(p.obj.j), 1)))
		indiv.pv <- rbind(indiv.pv, p.obj.j)
	}
	colnames(indiv.pv) <- paste0('D', 1:ncol(indiv.pv), '.p.value')
	rownames(indiv.pv) <- 1:nrow(indiv.pv)
	
	p.tab <- data.frame(indiv.pv, omni.p.value = omni.pv)
	rownames(p.tab) <- rownames(obj$aov.tab)[1:ncol(perm.mat)]
	
	
	list(p.tab = p.tab, aov.tab.list = aov.tab.list)
}



Rarefy <- function (otu.tab, depth = min(rowSums(otu.tab))){
	# Rarefaction function: downsample to equal depth
	#	
	# Args:
	#		otu.tab: OTU count table, row - n sample, column - q OTU
	#		depth: required sequencing depth 
	#
	# Returns:
	# 	otu.tab.rff: Rarefied OTU table
	#		discard: labels of discarded samples
	#
	otu.tab <- as.matrix(otu.tab)
	ind <- (rowSums(otu.tab) < depth)
	sam.discard <- rownames(otu.tab)[ind]
	otu.tab <- otu.tab[!ind, ]
	
	rarefy <- function(x, depth){
		y <- sample(rep(1:length(x), x), depth)
		y.tab <- table(y)
		z <- numeric(length(x))
		z[as.numeric(names(y.tab))] <- y.tab
		z
	}
	otu.tab.rff <- t(apply(otu.tab, 1, rarefy, depth))
	rownames(otu.tab.rff) <- rownames(otu.tab)
	colnames(otu.tab.rff) <- colnames(otu.tab)
	return(list(otu.tab.rff=otu.tab.rff, discard=sam.discard))
}


