  ## Generate families, plot alleles frequency distribution in the simulated population, plot number of genotype differences in the simulated population by comparing the individuals two by two - updated version
	
	# Data file
	datafile <- "input.csv"
	
	# Alleles and frequencies at each locus(individuals)
	matrix.indata <- scan(datafile, what = list(locus="",alleles=0,frequence=0),quiet=TRUE, sep=",",skip=1)
	
	# Number of individuals 
	nb.individus <- 1000000
	
	#Loci list and number
	locus.list <- c(matrix.indata$locus[1])
	for(i in 1:length(matrix.indata$locus))
	{
	if (matrix.indata$locus[i] != locus.list[length(locus.list)])
	{locus.list <- c(locus.list,matrix.indata$locus[i])}
	}
	nb.locus <- length(locus.list)
	
	# Matrix differences
	mat.diff = matrix(0, nrow=nb.individus, ncol=4, byrow=TRUE)
	colnames(mat.diff) <- c("Unrelated pair","Parent-offspring","Full sibs","Half sibs")
	
	#### N repetitions
	for (k in 1:nb.individus)
	{
	## Generate parents : P1, M, P2
	nb.parents = 3
	mat.parents = matrix(0, nrow=2*nb.parents, ncol=nb.locus, byrow=TRUE)
	rownames(mat.parents) <- c("P1_all-1","P1_all-2","M_all-1","M_all-2","P2_all-1","P2_all-2")
	colnames(mat.parents) <- locus.list
	
	for(i in 1:nb.locus){
	liste.alleles.locus <- c()
	freq.all.locus <- c()
	for(j in 1:length(matrix.indata$locus)){
	if(matrix.indata$locus[j] == locus.list[i]){
	liste.alleles.locus <- c(liste.alleles.locus, matrix.indata$alleles[j] )
	freq.all.locus <- c(freq.all.locus, matrix.indata$frequence[j] )
	}}
	mat.parents[,i] = sample(liste.alleles.locus, 2*nb.parents, prob=freq.all.locus, replace=TRUE) 
	}
	
	## Generate progeny : F1, F2, DF
	nb.progeny = 3
	mat.progeny = matrix(0, nrow=2*nb.progeny, ncol=nb.locus, byrow=TRUE)
	rownames(mat.progeny) <- c("F1_all-1","F1_all-2","F2_all-1","F2_all-2","DF_all-1","DF_all-2")
	colnames(mat.progeny) <- locus.list
	
	for(i in 1:nb.locus){
	mat.progeny[1,i] = sample(c(mat.parents[1,i],mat.parents[2,i]),1 , replace=TRUE ) # F1 allele pere
	mat.progeny[2,i] = sample(c(mat.parents[3,i],mat.parents[4,i]),1, replace=TRUE )# F1 allele mere
	mat.progeny[3,i] = sample(c(mat.parents[1,i],mat.parents[2,i]),1, replace=TRUE ) # F2 allele pere
	mat.progeny[4,i] = sample(c(mat.parents[3,i],mat.parents[4,i]),1, replace=TRUE ) # F2 allele mere 
	mat.progeny[5,i] = sample(c(mat.parents[5,i],mat.parents[6,i]),1, replace=TRUE ) # DF allele pere
	mat.progeny[6,i] = sample(c(mat.parents[3,i],mat.parents[4,i]),1, replace=TRUE ) # DF allele mere
	}
	
	## Calculate nb diff 
	nb.diff.ind.nonerelated = 0
	nb.diff.ind.pardesc = 0
	nb.diff.ind.pleinfrere = 0
	nb.diff.ind.demifrere = 0
	
	for(i in 1:nb.locus){
	
		# between P1 and P2 (none related)
	all.locP1 = c(mat.parents[1,i],mat.parents[2,i]) #allele list P1
	all.locP2 = c(mat.parents[5,i],mat.parents[6,i] )#allele list P2
	
	if(length(unique(c(all.locP1,all.locP2))) == 1){
	nb.diff.loc.nonerelated = 0
	}
	else{
	nb.diff.loc.nonerelated = 2 - length(intersect(all.locP1,all.locP2))
	}
	nb.diff.ind.nonerelated = nb.diff.ind.nonerelated + nb.diff.loc.nonerelated
	
		# between P1 and F1 (parent-offspring) 
	all.locF1 =  c(mat.progeny[1,i], mat.progeny[2,i]) #allele list F1
	
	if(length(unique(c(all.locP1,all.locF1))) == 1){
	nb.diff.loc.pardesc = 0
	}
	else{
	nb.diff.loc.pardesc = 2 - length(intersect(all.locP1,all.locF1))
	}
	nb.diff.ind.pardesc = nb.diff.ind.pardesc + nb.diff.loc.pardesc
	
		# between F1 and F2 (full siblings)
	all.locF2 =  c(mat.progeny[3,i], mat.progeny[4,i]) #allele list F2
	
	if(length(unique(c(all.locF1,all.locF2))) == 1){
	nb.diff.loc.pleinfrere = 0
	}
	else{
	nb.diff.loc.pleinfrere = 2 - length(intersect(all.locF1,all.locF2))
	}
	nb.diff.ind.pleinfrere = nb.diff.ind.pleinfrere + nb.diff.loc.pleinfrere
		
		# between F1 and HF (half siblings)
	all.locDF =  c(mat.progeny[5,i], mat.progeny[6,i]) #allele list HF
	
	if(length(unique(c(all.locF1,all.locDF))) == 1){
	nb.diff.loc.demifrere = 0
	}
	else{
	nb.diff.loc.demifrere = 2 - length(intersect(all.locF1,all.locDF))
	}
	nb.diff.ind.demifrere = nb.diff.ind.demifrere + nb.diff.loc.demifrere
	}
	
	mat.diff[k,] = c(nb.diff.ind.nonerelated, nb.diff.ind.pardesc, nb.diff.ind.pleinfrere, nb.diff.ind.demifrere)
	}
	
	## Plot the 4 distributions 
	tab.freq = matrix(0, nrow=4, ncol=2*nb.locus+1, byrow=TRUE)
	rownames(tab.freq) <- c("Unrelated pair","Parent-offspring","Full sibs","Half sibs")
	colnames(tab.freq) <- seq(0,2*nb.locus)
	
	for(i in 1:4)
	{
	for(j in 0:(2*nb.locus))
	{
	tab.freq[i,j+1] <- sum(mat.diff[,i]==j)/nb.individus
	}}
	write.csv(tab.freq,file="output.csv",row.names = TRUE)
	
	par(mfrow=c(2,2))
	barplot(tab.freq[1,], main = "Unrelated pair", xlab = 'Number of differences', ylab = 'Frequency',xlim = c(0,(2*nb.locus+2)),ylim = c(0,0.5)) 
	abline(h=0.05, col="blue")
	barplot(tab.freq[2,], main = "Parent-offspring", xlab = 'Number of differences', ylab = 'Frequency',xlim = c(0,(2*nb.locus+2)),ylim = c(0,0.5)) 
	abline(h=0.05, col="blue")
	barplot(tab.freq[3,], main = "Full sibs", xlab = 'Number of differences', ylab = 'Frequency',xlim = c(0,(2*nb.locus+2)),ylim = c(0,0.5)) 
	abline(h=0.05, col="blue")
	barplot(tab.freq[4,], main = "Half sibs", xlab = 'Number of differences', ylab = 'Frequency',xlim = c(0,(2*nb.locus+2)),ylim = c(0,0.5)) 
	abline(h=0.05, col="blue")
	
	
	## Calculus of heterozygozy
	he.distr <- c()
	for(i in 1:nb.locus) #loop on the loci
	{
	list.freq <- c()
	for(j in 1:length(matrix.indata$locus)){ #loop on the data file
	if(matrix.indata$locus[j] == locus.list[i]){
	list.freq <- c(list.freq, matrix.indata$frequence[j])
	}}
	som.freq.square <- sum(list.freq**2)
	he.distr <- c(he.distr, 1 - som.freq.square)
	}
	mat.he = matrix(he.distr, nrow=1, ncol=nb.locus, byrow=TRUE)
	rownames(mat.he) <- "he"
	colnames(mat.he) <- locus.list
	sorted.he.distr = sort(mat.he[1,], decreasing=TRUE)
	x11()
	barplot(sorted.he.distr, main = "Heterozygozy at each loci", xlab = 'Loci', ylab = 'Heterozygozy', names.arg=colnames(sorted.he.distr), axis.lty = 0, axisnames = TRUE, ylim = c(0,1))


