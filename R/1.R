set.seed(123)
# devtools::install_github("devproglab/BGVAR")
library(readxl)
library(tidyverse)
library(BGVAR)
library(glue)
options(scipen=9999)
path_raw <- '../../Science/Regional-Multipliers-RU/Data/РегДанные.xlsx'
path_data <- '../../Science/Regional-Multipliers-RU/Results/ModelData/'
path_result <- '../../Science/Regional-Multipliers-RU/Results/'
ncores <- 15
flag <- ''
saver <- function(object, folder, name) {
  saveRDS(object, paste0(folder, '/', name), compress = FALSE)
  gc()
}
codes <- read_excel(path_raw, sheet = 'КОДЫ') %>%
  select(region=name_official, abbrev=abbrev, OKATO_id) %>%
  arrange(abbrev)
# Preliminaries ----
### Weighing Matrix ----
W <- readRDS(glue('{path_data}W.RDS'))
### GRP Weights for Aggregation ----
aggW <- read_excel(path_raw, sheet='ФАКТОРЫ СВОД') %>%
  left_join(codes) %>%
  select(abbrev, w = GRPShareSample) %>%
  arrange(abbrev)
### Data ----
data <- readRDS(glue('{path_data}data.RDS'))
### Check Order  ----
all(colnames(W) == names(data[-80]))
all(rownames(W) == names(data[-80]))
all(colnames(W) == aggW[,1])
## Dominant Unit Settings ----
RU.weights <- list()
RU.weights$weights <- aggW %>% pull(w)
names(RU.weights$weights) <- aggW %>% pull(abbrev)
RU.weights$variables <- c('poil', 'expfed', 'e', 'revfed', 'i', 'pi', 'y')
RU.weights$exo <- c('poil', 'expfed', 'e', 'revfed', 'i')
OE.weights <- list(RU=RU.weights)
# Estimation ----
nvar <- 3
draws <- 200
burnin <- 0
thin = 2
plag <- c(1,1)
prior <- 'MN'
SV <- TRUE
mname <- glue("{nvar}v{paste0(plag,collapse='_')}p-{prior}-{ifelse(SV, 'SV', 'NSV')}-{draws}-{burnin}{flag}")
path_model <- glue('{path_result}{mname}/')
dir.create(path_model, showWarnings = FALSE)
if (file.exists(glue('{path_model}{mname}.RDS')) & !exists('model.1')) {
  model.1 <- readRDS(glue('{path_model}{mname}.RDS'))
} else if (!file.exists(glue('{path_model}{mname}.RDS'))) {
  model.1 <- bgvar(Data=data,
                   W = W,
                   draws = draws,
                   burnin = burnin,
                   plag = plag,
                   prior = prior,
                   hyperpara = list(prmean = ifelse(flag=='diff', 0, 1)),
                   SV = T,
                   thin = thin,
                   trend = T,
                   hold.out = 0,
                   eigen = 1.05,
                   expert = list(OE.weights=OE.weights,
                                 cores = ncores
                   )
  )
  saver(model.1, path_model, glue('{mname}.RDS'))
  msumm <- summary(model.1)
  saver(list(
    CD = msumm$CD$perc,
    res = msumm$res$p.res,
    crosscorr = msumm$cross.corr$res.res,
    crosscorrdata = msumm$cross.corr$dat.res
  ), path_model, glue('{mname}-summary.RDS'))
}
# IRF Calculation ----
QQ <- c(0.16, 0.5, 0.84)
## Federal Shocks ----
### Recursive Identification ----
ident <- 'chol'
shockinfo <- get_shockinfo(ident)
shockinfo$shock <- "RU.expfed"
shockinfo$scale <- 1
i3 <- irf(model.1, n.ahead=24, shockinfo=shockinfo, quantiles = QQ,
          expert = list(use_R = F, cores = ncores))
plot(i3, cumulative = F, resp = 'RU', quantiles = QQ)

# library(bigmemory)
# cores <- 10
# applyfun <- NULL
# if(is.null(applyfun)) {
#   applyfun <- if(is.null(cores)) {
#     lapply
#   } else {
#     if(.Platform$OS.type == "windows") {
#       cl_cores <- parallel::makeCluster(cores)
#       on.exit(parallel::stopCluster(cl_cores))
#       function(X, FUN, ...) parallel::parLapply(cl = cl_cores, X, FUN, ...)
#     } else {
#       function(X, FUN, ...) mcprogress::pmclapply(X, FUN, ..., mc.cores =
#                                                     cores)
#     }
#   }
# }
# irf_bigmat <- big.matrix(nrow=1000, ncol=1000,
#                          type="double", backingfile="IRF_data.bin",
#                          descriptorfile="IRF_data.desc")
# 
# imp.obj <- applyfun(1:1000,function(irep){
#   x <- bigmemory::attach.big.matrix('IRF_data.desc')
#   x[,irep] <- seq(1,1000)
# })
# o <- irf_bigmat[]
# 
# 
# library(filematrix)
# fm = fm.create(
#   filenamebase = "big_fm",
#   nrow = 1e5,
#   ncol = 1e5)
# cores <- 10
# applyfun <- NULL
# if(is.null(applyfun)) {
#   applyfun <- if(is.null(cores)) {
#     lapply
#   } else {
#     if(.Platform$OS.type == "windows") {
#       cl_cores <- parallel::makeCluster(cores)
#       on.exit(parallel::stopCluster(cl_cores))
#       function(X, FUN, ...) parallel::parLapply(cl = cl_cores, X, FUN, ...)
#     } else {
#       function(X, FUN, ...) mcprogress::pmclapply(X, FUN, ..., mc.cores =
#                                                     cores)
#     }
#   }
# }
# tic = proc.time()
# imp.obj <- applyfun(1:1000,function(irep){
#   fm = filematrix::fm.open(filenamebase = "big_fm", readonly = FALSE)
#   fm[,irep] = irep + 1:nrow(fm)
# })
# toc = proc.time()
# show(toc-tic)
# 
# fm[,4]
# 
# # tic = proc.time()
# # for( i in seq_len(ncol(fm)) ) {
# #   message(i, " of ", ncol(fm))
# #   fm[,i] = i + 1:nrow(fm)
# # }
# # toc = proc.time()
# # show(toc-tic)
# 
# # Cleanup
# 
# closeAndDeleteFiles(fm)
# 
