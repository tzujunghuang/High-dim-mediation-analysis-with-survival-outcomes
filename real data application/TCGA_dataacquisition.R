##### Configuration
# Make sure we have every package needed
package_list = c('survival', 'MASS', 'Matrix', 'mvtnorm', 'dplyr', 'stringr')

for (package in package_list) {
  require(package, character.only = TRUE); library(package, character.only = TRUE) }


#set working directory
#setwd('C:/Users/thuang2/maxMediSurv/R codes/real data application/')
source('data_management.R')

pheno <- read.csv('./clinical_754.csv', header = TRUE)

U <- pheno[, c('age_at_initial_pathologic_diagnosis', 'pathologic_stage', 
               'radiation_therapy', 'gender')] %>%
  rename(age = age_at_initial_pathologic_diagnosis,
         stage = pathologic_stage, radiation_idx = radiation_therapy) %>%
  mutate(stage = case_when(
           stage == "Stage I" ~ '0', 
           stage == "Stage IA" ~ '1',
           stage == "Stage IB" ~ '2', 
           stage == "Stage II" ~ '3',
           stage == "Stage IIA" ~ '4', 
           stage == "Stage IIB" ~ '5',
           stage == "Stage III" ~ '6', 
           stage == "Stage IIIA" ~ '7',
           stage == "Stage IIIB" ~ '8', 
           stage == "Stage IV" ~ '9',
           .default = as.character(stage)),
         stage = as.numeric(stage),
         radiation_idx = ifelse(radiation_idx == 'YES', 1, 0),
         gender = ifelse(gender == 'MALE', 1, 0)) %>%
  mutate(stage_f = factor(stage, levels = c(0:9), ordered = T),
         radiation_idx_f = factor(radiation_idx, levels = c(0,1), ordered = T),
         gender_f = factor(gender, levels = c(0,1), ordered = T)
         ) %>%
  mutate(stage_reduced = case_when(
           stage %in% c(0:2) ~ 0,
           stage %in% c(3:5) ~ 1,
           stage %in% c(6:9) ~ 2,
           .default = NA)) %>%
  mutate(stage_reduced_f = factor(stage_reduced, levels = c(0:2), ordered = T)) 
  # %>% select(age, stage_reduced_f, radiation_idx_f, gender_f)

print('Pheno data done!')


meth <- read.table('./TCGA.LUNG.sampleMap_HumanMethylation450', header = TRUE)
meth <- na.omit(meth)

p <- nrow(meth)
M <- meth[1:p, 2:ncol(meth)]
M <- apply(M, 2, as.numeric)
M <- as.data.frame(t(M))
colnames(M) <- as.vector(meth[1:p,1])

M_colSDs = colVars(M, sd_use = TRUE)
M = colStandardization(M, M_colSDs)

pheno$sampleID <- gsub('-','\\.',pheno$sampleID)
common_id <- intersect(pheno$sampleID, rownames(M))  # 754 common id, the same as pheno

M <- M[common_id,]    # sort to match the id
n <- nrow(M); p <- ncol(M)
print(c(n, p))

##### Showcase of partial data of mediators 
#save(M, file = 'Partial TCGA methylation.RData')


dat = list(X = log(pheno$OS), delta = pheno$Death, A = pheno$X.smoking..current.smoker_1.0., 
           B = M, U = U)

save(dat, file = 'TCGA lungcancer data.RData')


##### Temp
source('code_maxMediSurv.R')
X = dat$X; delta = dat$delta; A = dat$A; 

inverse_weight = Inverse_weight_func(input_data = dat, x0 = X, 
                                     obs = 1:n, quar_trunc = 0.9, 
                                     err_msg = 'Error in Est_Psi0_d Inverse_weight KM_table')
Y_Ghat = X*delta / inverse_weight; 
Y_Ghat[is.nan(Y_Ghat) | is.na(Y_Ghat)] = 0

dat0 = data.frame(Y_Ghat, A, dat$U)
t0 = proc.time()
pre.fit0 = lm(Y_Ghat ~ age + stage_reduced_f + 
                       radiation_idx_f + gender_f, data = dat0)
proc.time() - t0

U_mat = as.matrix(dat$U)

t1 = proc.time()
pre.fit0.alt = .lm.fit(U_mat, Y_Ghat)
proc.time() - t1

summary(pre.fit0)$coefficient

sig.confounders = do.call(rbind, 
    lapply(c('age', 'stage_reduced_f', 
             'radiation_idx_f', 'gender_f'), function(item){
    string.formula = paste0('Y_Ghat ~ ', item)
    pre.fit.marginal = lm(string.formula, data = dat0)
    
    if (summary(pre.fit.marginal)$coefficient[,4][2] < 0.1) {
      return( data.frame('confounder' = item, 
                         'p-val' = summary(pre.fit.marginal)$coefficient[,4][2]) )
    } }))
