
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


adonis3 <- function (formula, data = NULL, permutations = 999, method = "bray", 
		strata = NULL, contr.unordered = "contr.sum", contr.ordered = "contr.poly", 
		parallel = getOption("mc.cores"), ...)  {
	EPS <- sqrt(.Machine$double.eps)
	TOL <- 1e-07
	Terms <- terms(formula, data = data)
	lhs <- formula[[2]]
	lhs <- eval(lhs, data, parent.frame())
	formula[[2]] <- NULL
	rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
	op.c <- options()$contrasts
	options(contrasts = c(contr.unordered, contr.ordered))
	rhs <- model.matrix(formula, rhs.frame)
	options(contrasts = op.c)
	grps <- attr(rhs, "assign")
	qrhs <- qr(rhs)
	rhs <- rhs[, qrhs$pivot, drop = FALSE]
	rhs <- rhs[, 1:qrhs$rank, drop = FALSE]
	grps <- grps[qrhs$pivot][1:qrhs$rank]
	u.grps <- unique(grps)
	nterms <- length(u.grps) - 1
	if (nterms < 1) 
		stop("right-hand-side of formula has no usable terms")
	H.s <- lapply(2:length(u.grps), function(j) {
				Xj <- rhs[, grps %in% u.grps[1:j]]
				qrX <- qr(Xj, tol = TOL)
				Q <- qr.Q(qrX)
				tcrossprod(Q[, 1:qrX$rank])
			})
	if (inherits(lhs, "dist")) {
		if (any(lhs < -TOL)) 
			stop("dissimilarities must be non-negative")
		dmat <- as.matrix(lhs^2)
	}
	else if ((is.matrix(lhs) || is.data.frame(lhs)) && isSymmetric(unname(as.matrix(lhs)))) {
		dmat <- as.matrix(lhs^2)
		lhs <- as.dist(lhs)
	}
	else {
		dist.lhs <- as.matrix(vegdist(lhs, method = method, ...))
		dmat <- dist.lhs^2
	}
	n <- nrow(dmat)
	G <- -sweep(dmat, 1, rowMeans(dmat))/2
	SS.Exp.comb <- sapply(H.s, function(hat) sum(G * t(hat)))
	SS.Exp.each <- c(SS.Exp.comb - c(0, SS.Exp.comb[-nterms]))
	H.snterm <- H.s[[nterms]]
	tIH.snterm <- t(diag(n) - H.snterm)
	
	
	G.s <- list()
	G.s[[1]] <- G
	if (length(H.s) > 1)
		for (i in 2:(length(H.s))) 	G.s[[i]] <- (diag(n) - H.s[[i - 1]]) %*% G %*% (diag(n) - H.s[[i - 1]])
	
	
	if (length(H.s) > 1) 
		for (i in length(H.s):2) H.s[[i]] <- H.s[[i]] - H.s[[i - 1]]
	SS.Res <- sum(G * tIH.snterm)
	df.Exp <- sapply(u.grps[-1], function(i) sum(grps == i))
	df.Res <- n - qrhs$rank
	if (inherits(lhs, "dist")) {
		beta.sites <- qr.coef(qrhs, as.matrix(lhs))
		beta.spp <- NULL
	}
	else {
		beta.sites <- qr.coef(qrhs, dist.lhs)
		beta.spp <- qr.coef(qrhs, as.matrix(lhs))
	}
	colnames(beta.spp) <- colnames(lhs)
	colnames(beta.sites) <- rownames(lhs)
	F.Mod <- (SS.Exp.each/df.Exp)/(SS.Res/df.Res)
	f.test <- function(tH, G, df.Exp, df.Res, tIH.snterm) {
		(sum(G * tH)/df.Exp)/(sum(G * tIH.snterm)/df.Res)
	}
	SS.perms <- function(H, G, I) {
		c(SS.Exp.p = sum(G * t(H)), S.Res.p = sum(G * t(I - H)))
	}
	getPermuteMatrix <- getFromNamespace("getPermuteMatrix", "vegan")
	p <- getPermuteMatrix(permutations, n, strata = strata)
	permutations <- nrow(p)
	if (permutations) {
		tH.s <- lapply(H.s, t)
		if (is.null(parallel)) 
			parallel <- 1
		hasClus <- inherits(parallel, "cluster")
		isParal <- hasClus || parallel > 1
		isMulticore <- .Platform$OS.type == "unix" && !hasClus
		if (isParal && !isMulticore && !hasClus) {
			parallel <- makeCluster(parallel)
		}
		if (isParal) {
			if (isMulticore) {
				f.perms <- sapply(1:nterms, function(i) unlist(mclapply(1:permutations, 
											function(j) f.test(tH.s[[i]], G.s[[i]][p[j, ], p[j, 
																]], df.Exp[i], df.Res, tIH.snterm), mc.cores = parallel)))
			}
			else {
				f.perms <- sapply(1:nterms, function(i) parSapply(parallel, 
									1:permutations, function(j) f.test(tH.s[[i]], 
												G.s[[i]][p[j, ], p[j, ]], df.Exp[i], df.Res, tIH.snterm)))
			}
		}
		else {
			f.perms <- sapply(1:nterms, function(i) sapply(1:permutations, 
								function(j) f.test(tH.s[[i]], G.s[[i]][p[j, ], p[j, 
													]], df.Exp[i], df.Res, tIH.snterm)))
		}
		if (isParal && !isMulticore && !hasClus) 
			stopCluster(parallel)
		P <- (rowSums(t(f.perms) >= F.Mod - EPS) + 1)/(permutations + 
					1)
	}
	else {
		f.perms <- P <- rep(NA, nterms)
	}
	SumsOfSqs = c(SS.Exp.each, SS.Res, sum(SS.Exp.each) + SS.Res)
	tab <- data.frame(Df = c(df.Exp, df.Res, n - 1), SumsOfSqs = SumsOfSqs, 
			MeanSqs = c(SS.Exp.each/df.Exp, SS.Res/df.Res, NA), F.Model = c(F.Mod, 
					NA, NA), R2 = SumsOfSqs/SumsOfSqs[length(SumsOfSqs)], 
			P = c(P, NA, NA))
	rownames(tab) <- c(attr(attr(rhs.frame, "terms"), "term.labels")[u.grps], 
			"Residuals", "Total")
	colnames(tab)[ncol(tab)] <- "Pr(>F)"
	attr(tab, "heading") <- c("Terms added sequentially (first to last)\n")
	class(tab) <- c("anova", class(tab))
	out <- list(aov.tab = tab, call = match.call(), coefficients = beta.spp, 
			coef.sites = beta.sites, f.perms = f.perms, model.matrix = rhs, 
			terms = Terms)
	class(out) <- "adonis"
	out
}



PermanovaG <- function(formula, data = NULL, ...) {
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
	lhs <- eval(lhs, data, parent.frame())
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
		
		obj <- adonis(as.formula(paste("Y", "~", rhs)), data, ...)
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



PermanovaG2 <- function(formula, data = NULL, ...) {
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
	lhs <- eval(lhs, data, parent.frame())
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
		
		obj <- adonis3(as.formula(paste("Y", "~", rhs)), data, ...)
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

calculateK <- function (G, H, df2) {
	n <- nrow(G)
	dG <- diag(G)
	
	mu1 <- sum(dG) / (n - df2)
	mu2 <- (sum(G^2) - sum(dG^2)) / ((n - df2)^2 + sum(H^4) - 2 * sum(diag(H^2)))
	
	k <- mu1^2 / mu2
	list(K = k, mu1 = mu1, mu2 = mu2)
}



dmanova <- function (formula, data = NULL, positify = FALSE,
		contr.unordered = "contr.sum", contr.ordered = "contr.poly", 
		returnG = FALSE) {
	Terms <- terms(formula, data = data)
	lhs <- formula[[2]]
	lhs <- eval(lhs, data, parent.frame())
	
	formula[[2]] <- NULL
	rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
	op.c <- options()$contrasts
	options(contrasts = c(contr.unordered, contr.ordered))
	rhs <- model.matrix(formula, rhs.frame)
	options(contrasts = op.c)
	grps <- attr(rhs, "assign")
	
	# The first group includes the intercept, nterms count the intercept
	u.grps <- unique(grps)
	nterms <- length(u.grps)
	
	Z <- rhs[, grps %in% u.grps[1:(nterms-1)], drop=F]
	XZ <- rhs[, grps %in% u.grps[1:(nterms)], drop=F]
	
	n <- nrow(XZ)
	if (class(lhs) == 'dist') {
		D <- as.matrix(lhs)
	} else {
		D <- lhs
	}
	
	D <- -D^2 / 2
	
	G <- mean(D) + D - rowMeans(D) - matrix(rep(1, n), ncol=1) %*% colMeans(D)
	
	if (positify) {
		G <- as.matrix(nearPD(G)$mat)
	}
	
	XZi <- solve(t(XZ) %*% XZ)
	Zi <- solve(t(Z) %*% (Z))
	HZ <- Z %*% Zi %*% t(Z)
	HXZ <- XZ %*% XZi %*% t(XZ)
	HX <- HXZ - HZ
	HIXZ <- diag(n) - HXZ
	HIX <- diag(n) - HX
	
	
	df1 <- ncol(XZ) - ncol(Z)
	df2 <- n - ncol(XZ)
	
	MSS <- sum(G * HX)
	RSS <- sum(G * HIXZ)
	TSS <- sum(diag(G))
	
	
	f.stat <- (MSS / df1) / (RSS / df2)
	
	GXZ <- G %*% XZ
	XZXZi <- XZ %*% XZi
	GXZtXZXZi <- GXZ %*% t(XZXZi)
	
	G.tilde <- G + XZXZi %*% (t(GXZ) %*% XZ) %*% t(XZXZi) - GXZtXZXZi - t(GXZtXZXZi)
	
		
	obj <- calculateK(G.tilde, HIXZ, n - df2)
	K <- obj$K
	
	
	p.value <- pchisq(f.stat * K * df1, df = K * df1, lower.tail = F)
		


	SumsOfSqs <- c(MSS, RSS, TSS)
	
	tab <- data.frame(Df = c(df1, df2, n - 1), SumsOfSqs = SumsOfSqs, 
			MeanSqs = c(MSS/df1, RSS/df2, NA), F.Model = c(f.stat, 
					NA, NA), R2 = c(MSS/TSS, NA, NA), 
			"Pr(>F)" = c(p.value, NA, NA))
	
	rownames(tab) <- c("Model(Adjusted)", "Residuals", "Total")
	colnames(tab)[ncol(tab)] <- c("Pr(>F)")
	class(tab) <- c("anova", class(tab))
	attr(tab, "heading") <- c("F stat and P value of the last term is adjusted by preceding terms!\n")
	
	if (returnG == TRUE) {
		out <- list(aov.tab = tab,  df = K, G = G,  call = match.call())
	} else {
		out <- list(aov.tab = tab, df = K, call = match.call())
	}
	
	out
}