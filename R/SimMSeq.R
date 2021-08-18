

#' @title estimiate the DM parameters
#' @param alpha shape factor
#' @return a matrix represents the true proportion of the ref.otu.tab, row is otu, column is sample
#' @import stats
#' @rdname rdirichlet.m
#' @export

rdirichlet.m <- function (alpha) {
	Gam <- matrix(rgamma(length(alpha), shape = alpha), nrow(alpha), ncol(alpha))
	t(t(Gam) / colSums(Gam))
}



#' @title estimiate the true proportion of ref.otu.tab
#' @param ref.otu.tab shape factor
#' @return A list with the elements
#' \item{mu}{absolute abundance before applying size factor}
#' \item{ref.otu.tab}{absolute abundance renormalized by size factor}
#' @rdname EstPara
#' @import dirmult
#' @import stats
#' @export

EstPara <- function (ref.otu.tab) {
	
	if (is.null(rownames(ref.otu.tab))) {
		rownames(ref.otu.tab) <- paste0('OTU', 1 : nrow(ref.otu.tab))
	} # otu * sample
	samplenames = colnames(ref.otu.tab)
	taxnames = rownames(ref.otu.tab)
	
	dirmult.paras <- dirmult::dirmult(t(ref.otu.tab))
	
	gamma = dirmult.paras$gamma
	names(gamma) = names(dirmult.paras$pi)
	
	# Add pseduo count(each OTU add gamma estimated from dirmult)
	ref.otu.tab = sapply(1:ncol(ref.otu.tab), function (i) gamma + ref.otu.tab[,i]) # C_ij otu * sample
	
	# back to dirchlet, calculate the true proportion
	ref.otu.tab.p <- rdirichlet.m(ref.otu.tab) # P_ij nOTU*nSam
	colnames(ref.otu.tab.p) = samplenames
	rownames(ref.otu.tab.p) = taxnames
	
	# order OTUs by mean OTU proportion, for later selection
	ord = order(rowMeans(ref.otu.tab.p), decreasing = TRUE)
	ref.otu.tab.p =  ref.otu.tab.p[ord,]
	
	# apply size factor
	Si = exp(rnorm(ncol(ref.otu.tab.p)))
	ref.otu.tab0 = t(t(ref.otu.tab.p)*Si)
	colnames(ref.otu.tab0) = colnames(ref.otu.tab.p)
	return(list(mu = ref.otu.tab.p, ref.otu.tab = ref.otu.tab0))
}





#' @title Simulate microbiome data for differential abundance analysis
#' @param ref.otu.tab taxa in rows, and samples in columns.
#' @param nSam an integer: indicating the number of samples you want to simulate.
#' @param nOTU an integer: indicating the number of taxa you want to simulate.
#' @param diff.otu.pct a number between 0-1: indicating the proportion of differential taxa, default is 0.
#' @param diff.otu.direct a character: consists "balanced" and "unbalanced". "balanced" is default, means differential taxa increase and decrease, "unbalanced" means differetial taxa only increase or only decrease.
#' @param diff.otu.mode a character: consists "rare", "mix", "abundant". "rare" means differential taxa come from the least abundant 0.25 taxa. "abundant", is default, means the most abundant 0.25 taxa. "mix" means 0.5 of the differential taxa come from "abundnat" and "rare", respectively.
#' @param covariate.type a character:"binary", "continuous". "binary" is default, means covariate is binary variable, i.e., case vs control. "continuous" means the covariate is continuous variable, i.e., age.
#' @param covariate.eff.mean a number: indicating the effect size, default is 0.
#' @param grp.ratio a number between 0-1: indicating the ratio of case and control. default is 1, means 50 percentage samples are case and  50 percentage  are control. formula is grp.ratio / (1 + grp.ratio)).
#' @param covariate.eff.sd a number: indicating the variation of the taxa. default is 0, default is 0.
#' @param confounder.type a character: consists "none", "binary", "continuous", "both". default is "none", means no confounding variable. binary" means confounder is 2-group variable, i.e., sex. "continuous" means a continuous variable, i.e., BMI. "both" means two confouders, including binary and continuous variable.
#' @param conf.cov.cor a number between 0-1: indicating the correlation between covariate and confounders, default is 0.6.
#' @param conf.diff.otu.pct a number between 0-1: indicating the proportion of taxa are influenced by confounders, default is 0.05.
#' @param conf.nondiff.otu.pct a number between 0-1: indicating the proportion of taxa are not influenced by confounders, default is 0.1.
#' @param confounder.eff.mean a number: effect size of the confounding effect, default is 0.
#' @param confounder.eff.sd a number: variation of the taxa that have confounding effect. default is 0.
#' @param depth.mu an integer: define the mean sequencing depth to be simulated, default is 10000.
#' @param depth.theta theta value in defining the variance 'mu + mu^2/theta' of the negative binomial distribution, default is 5.
#' @param depth.conf.factor a number: depth confounding factor, 0 means no depth confounding, default is 0
#' @return Return a list with the elements:
#' \item{otu.tab.sim}{simulated otu table}
#' \item{covariate}{simulated covariate}
#' \item{diff.otu.ind}{index of the differential taxa}
#' \item{otu.names}{the otu names of the simulated otu table}
#' \item{conf.otu.ind}{taxa with confounding effect}
#' @rdname SimulateSeq
#' @import dirmult
#' @import MASS
#' @import stats
#' @export

## Little different from the general simulation function in cluster, since we want to make the nOTU and nSam is the same as the nrow/ncol of the input ref.otu.tab
SimulateMSeq <- function (
		ref.otu.tab, nSam = 100, nOTU = 500,
		diff.otu.pct = 0.1, diff.otu.direct = c("balanced", "unbalanced"), diff.otu.mode = c("abundant", "rare", "mix"),
		covariate.type = c("binary", "continuous"), grp.ratio = 1, covariate.eff.mean = 1, covariate.eff.sd = 0,
		confounder.type = c("none", "binary", "continuous", "both"), conf.cov.cor = 0.6, 
		conf.diff.otu.pct = 0, conf.nondiff.otu.pct = 0.1, confounder.eff.mean = 0, confounder.eff.sd = 0,
		error.sd = 0, depth.mu = 10000, depth.theta = 5,  depth.conf.factor = 0
){
	diff.otu.direct <- match.arg(diff.otu.direct)
	diff.otu.mode <- match.arg(diff.otu.mode)
	covariate.type <- match.arg(covariate.type)
	confounder.type <- match.arg(confounder.type)
	ref.otu.tab0 <- ref.otu.tab
	model.paras <- EstPara(ref.otu.tab = ref.otu.tab) # full ref.otu.tab + estimated parameters
	
	sample.names <- colnames(model.paras$ref.otu.tab)
	ref.otu.tab <- model.paras$ref.otu.tab[(1 : (nOTU)),] # otu * sample
	dim(ref.otu.tab)
	idx.otu <- rownames(ref.otu.tab)
	idx.sample <- sample(sample.names, nSam)
	idx.nonsample <- colnames(ref.otu.tab)[!(colnames(ref.otu.tab) %in% idx.sample)]
	
	ref.otu.tab = ref.otu.tab[,idx.sample] # otu * sample, final selected OTU tab by User
	ref.otu.tab.unselect = ref.otu.tab0[c(1 : (nOTU)),][,idx.nonsample] # For overfitting test
	
	
	if (confounder.type == "none")  {
		confounder.type <- "continuous"
		confounder.eff.mean <- 0
		confounder.eff.sd <- 0
		Z <- NULL
	}
	if (confounder.type == "continuous") Z <- cbind(rnorm(nSam))
	if (confounder.type == "binary") Z <- cbind(c(rep(0, nSam %/% 2), rep(1, nSam - nSam %/% 2)))
	if (confounder.type == "both") Z <- cbind(rnorm(nSam), c(rep(0, nSam %/% 2), rep(1, nSam - nSam %/% 2)))
	
	rho <- sqrt(conf.cov.cor ^ 2 / (1 - conf.cov.cor ^ 2))
	
	if (covariate.type == 'continuous') {
		X <- rho * scale(scale(Z) %*% rep(1, ncol(Z)))  + rnorm(nSam)
	}
	
	if (covariate.type == "binary") {
		X <- rho * scale(scale(Z) %*% rep(1, ncol(Z)))  + rnorm(nSam)
		X <- cbind(ifelse(X <= quantile(X, grp.ratio / (1 + grp.ratio)), 0, 1))
	}
	rownames(X) <- colnames(ref.otu.tab)
	
	covariate.eff.mean1 = covariate.eff.mean
	covariate.eff.mean2 = covariate.eff.mean
	
	if (diff.otu.direct == "balanced") {
		if (diff.otu.mode == 'abundant'){
			eta.diff <- sample(c(rnorm(floor(nOTU / 2), mean = -covariate.eff.mean2, sd = covariate.eff.sd),
							rnorm(nOTU - floor(nOTU / 2), mean = covariate.eff.mean2, sd = covariate.eff.sd)))  %*% t(scale(X)) # multiple means * xb, eta = xb
		}else if(diff.otu.mode == 'rare'){
			eta.diff <- sample(c(rnorm(floor(nOTU / 2), mean = -covariate.eff.mean2, sd = covariate.eff.sd),
							rnorm(nOTU - floor(nOTU / 2), mean = covariate.eff.mean2, sd = covariate.eff.sd)))  %*% t(scale(X))
		}else{
			eta.diff <- c(sample(c(rnorm(floor(nOTU / 4), mean = -covariate.eff.mean1, sd = covariate.eff.sd),
									rnorm(floor(nOTU / 2) - floor(nOTU / 4), mean = covariate.eff.mean1, sd = covariate.eff.sd))),
					sample(c(rnorm(floor((nOTU -floor(nOTU / 2))/2), mean = -covariate.eff.mean2, sd = covariate.eff.sd),
									rnorm(nOTU -floor(nOTU / 2) - floor((nOTU -floor(nOTU / 2))/2), mean = covariate.eff.mean2, sd = covariate.eff.sd)))) %*% t(scale(X))
		}
	}
	
	
	if (diff.otu.direct == "unbalanced") {
		if (diff.otu.mode == 'abundant'){
			eta.diff <- rnorm(nOTU, mean = covariate.eff.mean2, sd = covariate.eff.sd) %*% t(scale(X)) # otu in 50% samples +, 50% sample -
		}else if(diff.otu.mode == 'rare'){
			eta.diff <- rnorm(nOTU, mean = covariate.eff.mean2, sd = covariate.eff.sd) %*% t(scale(X))
		}else{
			eta.diff <- c(sample(c(rnorm(floor(nOTU / 2), mean = covariate.eff.mean1, sd = covariate.eff.sd))),
					sample(c(rnorm(nOTU - floor(nOTU / 2), mean = covariate.eff.mean2, sd = covariate.eff.sd)))) %*% t(scale(X))
		}
	}
	
	eta.conf <-  sample(c(rnorm(floor(nOTU / 2), mean = -confounder.eff.mean, sd = confounder.eff.sd),
					rnorm(nOTU - floor(nOTU / 2), mean = confounder.eff.mean, sd = confounder.eff.sd)))  %*% t(scale(scale(Z) %*% rep(1, ncol(Z))))
	
	otu.ord <- 1 : (nOTU)
	diff.otu.ind <- NULL
	diff.otu.num <- round(diff.otu.pct * nOTU)
	
	if (diff.otu.mode == 'mix') diff.otu.ind <- c(diff.otu.ind, sample(otu.ord, diff.otu.num))
	if (diff.otu.mode == 'abundant') diff.otu.ind <- c(diff.otu.ind, sample(otu.ord[1 : round(length(otu.ord) / 4)], diff.otu.num))
	if (diff.otu.mode == 'rare') diff.otu.ind <- c(diff.otu.ind, sample(otu.ord[round(3 * length(otu.ord) / 4) : length(otu.ord)], diff.otu.num))
	
#	if(diff.otu.pct > 0){
#		conf.otu.ind <- c(sample(diff.otu.ind, round(nOTU * conf.diff.otu.pct)), sample(setdiff(1 : (nOTU), diff.otu.ind), round(conf.nondiff.otu.pct * nOTU)))
#		eta.diff[setdiff(1 : (nOTU), diff.otu.ind), ] <- 0 # some OTUs will be 0
#		eta.conf[setdiff(1 : (nOTU), conf.otu.ind), ] <- 0
#	}else{
#		eta.diff <- eta.diff
#		eta.conf <- eta.conf
#		conf.otu.ind <- NULL
#	}
	
    if (length(diff.otu.ind) >= round(nOTU * conf.diff.otu.pct)) {
		conf.otu.ind1 <- sample(diff.otu.ind, round(nOTU * conf.diff.otu.pct))
	} else {
		conf.otu.ind1 <- diff.otu.ind
	}
    conf.otu.ind <- c(conf.otu.ind1, sample(setdiff(1 : (nOTU), diff.otu.ind), round(conf.nondiff.otu.pct * nOTU)))
	
    eta.diff[setdiff(1 : (nOTU), diff.otu.ind), ] <- 0 # some OTUs will be 0
    eta.conf[setdiff(1 : (nOTU), conf.otu.ind), ] <- 0
	eta.error <- matrix(rnorm(nOTU * nSam, 0, error.sd), nOTU, nSam)

	
	eta.exp <- exp(t(eta.diff + eta.conf + eta.error))
	eta.exp <- eta.exp * t(ref.otu.tab)
	
	
	ref.otu.tab.prop <- eta.exp / rowSums(eta.exp)
	ref.otu.tab.prop <- t(ref.otu.tab.prop) # otu * sample
	
	nSeq <- rnegbin(nSam, mu = depth.mu * exp(scale(X) * depth.conf.factor), theta = depth.theta)

	otu.tab.sim <- sapply(1:ncol(ref.otu.tab.prop), function (i) rmultinom(1, nSeq[i], ref.otu.tab.prop[,i]))
	
	colnames(otu.tab.sim) <- rownames(eta.exp)
	rownames(otu.tab.sim) <- rownames(ref.otu.tab)
	
	
	diff.otu.ind = (1 : nOTU) %in% diff.otu.ind
	conf.otu.ind = (1 : nOTU) %in% conf.otu.ind
	
	return(list(otu.tab.sim = otu.tab.sim, covariate = X, confounder = Z, diff.otu.ind = diff.otu.ind, otu.names = idx.otu, conf.otu.ind = conf.otu.ind))
}
