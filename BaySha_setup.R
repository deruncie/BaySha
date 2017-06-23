library(data.table)

# Load the haplotye and map data, and then use R/qtl to impute missing genotypes
library(qtl)
haplotypes = do.call(rbind,lapply(1:5,function(x) read.delim(sprintf('SFP578map_At%s_genotypes.txt',x),skip=1,check.names = F,stringsAsFactors = F)))
rownames_haplotypes = haplotypes[,1]
haplotypes = haplotypes[,-1]
rownames(haplotypes) = rownames_haplotypes

map = do.call(rbind,lapply(1:5,function(x) data.frame(Chr=x,read.delim(sprintf('SFP578map_At%s_jm.txt',x),h=F,stringsAsFactors = F))))[,c(2,1,3)]
colnames(map) = c('Marker','Chr','cM')
write.csv(map,'Joint_map.csv')

stopifnot(all(rownames(haplotypes) %in% map[,1]))
haplotypes = haplotypes[map[,1],]

## build the cross file for R/qtl
cross_data = cbind(rbind(c('ID','',''),map[,1:3]),rbind(1:ncol(haplotypes),haplotypes))
colnames(cross_data) = NULL
write.csv(cross_data,file = 'RQtl_cross.csv',row.names = F)

# using R/qtl, read in the cross file, use fill.geno to impute
cross = read.cross('csvr','.','RQtl_cross.csv',genotypes = c('A','B'),na.strings = c('-'),crosstype = 'bc')
cross = convert2riself(cross)
cross <- calc.genoprob(cross, step=0)
cross = fill.geno(cross,method = 'argmax')

# collect genotypes, create genotype matrix
geno = do.call(rbind,lapply(cross$geno,function(x) t(x$data)))
rownames(geno) = map[,1]
colnames(geno) = colnames(haplotypes)

write.csv(geno,file = 'Rqtl_imputed_genotypes.csv')
geno = read.csv('Rqtl_imputed_genotypes.csv',check.names = F,row.names = 1)-1
geno = t(geno)


# load expression data and map lines to genotypes
sample_info = read.delim('E-TABM-126.sdrf.txt',stringsAsFactors=F)
expression = fread('E-TABM-126-processed-data-1355696613.txt',skip=1)
expression = as.data.frame(expression)
expression_header = fread('E-TABM-126-processed-data-1355696613.txt',nrows=1,h=F)
colnames(expression) = expression_header[1,]
colnames(expression) = sub('.CEL','',colnames(expression))

AffyID = expression[,1]
AffyID = sub('Affymetrix:CompositeSequence:ATH1-121501:','',AffyID)
affy_info = read.delim('A-AFFY-2.adf.txt',skip=14,stringsAsFactors = F)
affy_info = affy_info[match(AffyID,affy_info[,1]),]

# gene cis-QTLs
library(ath1121501.db)

# get probe info
# GeneID
mapped_probes <- mappedkeys(ath1121501ACCNUM)
affy_info$InDatabase = affy_info$Composite.Element.Name %in% mapped_probes
affy_info$Gene[affy_info$InDatabase] = do.call(c,as.list(ath1121501ACCNUM[affy_info$Composite.Element.Name[affy_info$InDatabase]]))

# Gene location
mapped_probes <- mappedkeys(ath1121501CHRLOC)
affy_info$InDatabase = affy_info$Composite.Element.Name %in% mapped_probes
locs = as.list(ath1121501CHRLOC[affy_info$Composite.Element.Name[affy_info$InDatabase]])
Chrs = sapply(locs,function(x) names(x)[1])
Poss = sapply(locs,function(x) abs(x[1]))
affy_info$Chr[affy_info$InDatabase] = Chrs
affy_info$Pos[affy_info$InDatabase] = Poss

# map chromosome positions
mapped_genes = mappedkeys(org.At.tairCHRLOC)
map$Gene = toupper(sapply(map[,1],function(x) strsplit(x,'-')[[1]][1]))
map$Gene_in_TAIR = map$Gene %in% mapped_genes
locs = as.list(org.At.tairCHRLOC[map$Gene[map$Gene_in_TAIR]])
Chrs = sapply(locs,function(x) names(x)[1])
Poss = sapply(locs,function(x) abs(x[1]))
map$Chr2[map$Gene_in_TAIR] = Chrs  # to check for consistency
map$Pos[map$Gene_in_TAIR] = Poss

# confirm that cM and Pos are related
with(subset(map,Chr==2),plot(cM~Pos))
# build spline models to predict cM from Pos
cM_mods = lapply(1:5,function(x) with(subset(map,Chr==x),approxfun(Pos,cM,rule=2)))
# predict cM from Pos for probes
for(x in 1:5){
  i = affy_info$Chr == x & !is.na(affy_info$Chr)
  affy_info$cM[i] = cM_mods[[x]](affy_info$Pos[i])
}

write.csv(affy_info,'extended_probe_annotations.csv')


probes_cross = data.frame(affy_info[affy_info$InDatabase,c('Composite.Element.Name','Chr','cM')],matrix(NA,sum(affy_info$InDatabase),ncol(haplotypes)))
colnames(probes_cross) = NULL
write.csv(probes_cross,file = 'RQtl_probes_cross.csv',row.names = F)
system('cat RQtl_cross.csv RQtl_probes_cross.csv > RQtl_joint_probes_cross.csv')


probe_imputation = read.cross('csvr','.','RQtl_joint_probes_cross.csv',genotypes = c('A','B'),na.strings = c('-','NA'),crosstype = 'bc')
probe_imputation = convert2riself(probe_imputation)
probe_imputation = fill.geno(probe_imputation,method = 'argmax')

probe_X = do.call(rbind,lapply(probe_imputation$geno,function(x) t(x$data))) - 1
probe_X = probe_X[rownames(probe_X) %in% affy_info$Composite.Element.Name,]
colnames(probe_X) = colnames(haplotypes)
write.csv(probe_X,file = 'Rqtl_imputed_cisGenotypes.csv',row.names = T)

