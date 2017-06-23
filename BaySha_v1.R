library(MCMCpack)
library(BSFG)
library(reshape2)
library(ggplot2)
library(data.table)


geno = read.csv('Rqtl_imputed_genotypes.csv',check.names = F,row.names = 1)-1
geno = t(geno)

# load expression data and map lines to genotypes
sample_info = read.delim('E-TABM-126.sdrf.txt',stringsAsFactors=F)
expression = fread('E-TABM-126-processed-data-1355696613.txt',skip=1)
expression = as.data.frame(expression,stringsAsFactors=F)
expression_header = fread('E-TABM-126-processed-data-1355696613.txt',nrows=1,h=F)
colnames(expression) = expression_header[1,]
colnames(expression) = sub('.CEL','',colnames(expression))

AffyID = expression[,1]
AffyID = sub('Affymetrix:CompositeSequence:ATH1-121501:','',AffyID)
expression = t(expression[,-1])
log_expression = log(expression)
colnames(log_expression) = AffyID

lines = read.csv('Lines.csv')

data = data.frame(Sample = rownames(log_expression))
sample_info = sample_info[match(data$Sample,sample_info$Hybridization.Name),]
data = cbind(data,Line = sample_info$Characteristics..StrainOrLine.)
data$ID = as.character(lines$Loudet.RIL...or.accession[match(data$Line,lines$ABRC.stock..)])


# factor QTL data
X_F = geno[data$ID,]

# gene cis-QTLs
probe_X = t(read.csv('Rqtl_imputed_cisGenotypes.csv',row.names = 1,check.names=F))
probe_X = probe_X[data$ID,]


probes = which(colnames(log_expression) %in% colnames(probe_X))[1:100]

Y = log_expression[,probes]
cis_genotypes = lapply(colnames(Y),function(x) probe_X[,x,drop=F])



# Set up BSFG
# initialize priors
run_parameters = BSFG_control(
  sampler = 'fast_BSFG',
  # sampler = 'general_BSFG',
  scale_Y = FALSE,
  simulation = FALSE,
  h2_divisions = 20,
  h2_step_size = NULL,
  burn = 00,
  k_init = 20
)

priors = list(
  # fixed_var = list(V = 5e5,   nu = 2.001),
  fixed_var = list(V = 1,     nu = 3),
  QTL_resid_var = list(V = 1/1000,     nu = 3),
  QTL_factors_var = list(V = 1/1000,     nu = 3),
  tot_Y_var = list(V = 0.5,   nu = 3),
  tot_F_var = list(V = 18/20, nu = 20),
  delta_1   = list(shape = 2.1,  rate = 1/20),
  delta_2   = list(shape = 3, rate = 1),
  Lambda_df = 3,
  B_df      = 3,
  B_F_df    = 3,
  # h2_priors_resids_fun = function(h2s,n) pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
  # h2_priors_factors_fun = function(h2s,n) ifelse(h2s == 0,n,n/(n-1))
  h2_priors_resids = 1,
  h2_priors_factors = 1
)


BSFG_state = BSFG_init(Y, model=~1+(1|ID), data, #factor_model_fixed = ~1,
                                  QTL_factors = X_F,
                                  priors=priors,run_parameters=run_parameters,
                                  cis_genotypes = cis_genotypes
                                  # K_mats = list(animal = K)
                          )

BSFG_state = reorder_factors(BSFG_state)
BSFG_state = clear_Posterior(BSFG_state)
n_samples = 100;
for(i  in 1:10) {
  print(sprintf('Run %d',i))
  BSFG_state = sample_BSFG(BSFG_state,n_samples,ncores=1)
  if(BSFG_state$current_state$nrun < BSFG_state$run_parameters$burn) {
    BSFG_state = reorder_factors(BSFG_state)
  }
  BSFG_state = save_posterior_chunk(BSFG_state)
  print(BSFG_state)
  plot(BSFG_state)
}

BSFG_state$Posterior = reload_Posterior(BSFG_state)


cross = read.cross('csvr','.','RQtl_cross.csv',genotypes = c('A','B'),na.strings = c('-'),crosstype = 'bc')
cross = convert2riself(cross)
cross = calc.genoprob(cross)

U_F = get_posterior_mean(BSFG_state,U_F,samples = 149:150)
rownames(U_F) = dimnames(BSFG_state$current_state$U_F)[[1]]
cross$pheno = U_F[rownames(geno),]
scanone_results = scanone(cross, phe = 1:ncol(U_F),method="hk")
sapply(1:(ncol(U_F)/3),function(x) plot(scanone_results,lodcolumn=1:3+(x-1)*3))

B = get_posterior_mean(BSFG_state,B_F %*% t(Lambda))
image(Matrix(B),aspect = 1)

