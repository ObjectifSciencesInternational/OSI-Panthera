## Generate individuals, plot alleles frequency distribution in the simulated population, plot number of genotype differences in the simulated population by comparing the individuals two by two - updated version

# Data file
datafile <- "input.csv"

# Alleles and frequencies at each locus(individuals)
matrix.indata <- scan(datafile, what = list(locus="",alleles=0,frequence=0),quiet=TRUE, sep=",",skip=1)

# Number of individuals 
nb.individus <- 1000

#Loci list and number
locus.list <- c(matrix.indata$locus[1])
for(i in 1:length(matrix.indata$locus))
{
if (matrix.indata$locus[i] != locus.list[length(locus.list)])
{locus.list <- c(locus.list,matrix.indata$locus[i])}
}
nb.locus <- length(locus.list)


# Generate matrix of individuals
matrix.all.ind = matrix(0, nrow=2*nb.individus, ncol=nb.locus)

savind <- c()
for (i in seq(from=1, to=(2*nb.individus-2),by=2)) # loop on the first sample
{
for (j in seq(from=i+2, to=(2*nb.individus),by=2)) # loop on the second sample
{
nb.diff.ind = 0
for (k in 1:nb.locus) #loop on all loci
{
liste.alleles.locus <- c()
freq.all.locus <- c()
for(l in 1:length(matrix.indata$locus)){
if(matrix.indata$locus[l] == locus.list[k]){
liste.alleles.locus <- c(liste.alleles.locus, matrix.indata$alleles[l] )
freq.all.locus <- c(freq.all.locus, matrix.indata$frequence[l] )
}}
matrix.all.ind[,k] = sample(liste.alleles.locus, 2*nb.individus, prob=freq.all.locus, replace=TRUE) 
#Counting of differences
all.loc1 = c(matrix.all.ind[i,k],matrix.all.ind[i+1,k]) #Allele list ind1
all.loc2 = c(matrix.all.ind[j,k], matrix.all.ind[j+1,k]) #Allele list ind2
if(length(unique(c(all.loc1, all.loc2))) == 1){
nb.diff.loc = 0
}
else{
nb.diff.loc = 2 - length(intersect(all.loc1,all.loc2))
}
nb.diff.ind = nb.diff.ind + nb.diff.loc
}
savind <-c(savind,nb.diff.ind)
}}
print(length(savind))

#Allele of differences
par(mfrow=c(3,3))
for(i in 1:nb.locus){
liste.alleles = sort(unique(matrix.all.ind[,i]))
count.alleles = c()
for(j in liste.alleles){
count.alleles = c(count.alleles, length(which(matrix.all.ind[,i] == j)))
}
plot(count.alleles, main = locus.list[i], xlab = 'Alleles', ylab = 'Number of occurrences', xaxt = "n",) 
axis(1, at=1:length(liste.alleles), labels=liste.alleles)
}

x11()
#print(savind)
hist(savind, main = "Number of differences between the genotypes for a simulated population of 1000 individuals", xlab = 'Number of differences', ylab = 'Number of occurrences') 

# Calculus of heterozygozy
he.distr <- c()
for(i in 1:nb.locus) #loop on the loci
{
som.freq.square <- 0
for(j in 1:length(matrix.indata$locus)){ #loop on the data file
if(matrix.indata$locus[j] == locus.list[i]){
som.freq.square = som.freq.square + matrix.indata$frequence[j]*matrix.indata$frequence[j]
}}
he.distr <- c(he.distr, 1 - som.freq.square)
}
mat.he = matrix(he.distr, nrow=1, ncol=nb.locus, byrow=TRUE)
rownames(mat.he) <- "he"
colnames(mat.he) <- locus.list
sorted.he.distr = sort(mat.he[1,], decreasing=TRUE)
print(sorted.he.distr)
x11()
barplot(sorted.he.distr, main = "Heterozygozy for each loci", xlab = 'Loci', ylab = 'Heterozygozy', names.arg=colnames(sorted.he.distr), axis.lty = 0, axisnames = TRUE, ylim = c(0,1))
