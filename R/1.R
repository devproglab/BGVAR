# set.seed(123)
# # devtools::install_github("devproglab/BGVAR", INSTALL_opts = c("--with-keep.source"), force = T)
# # devtools::install_github("devproglab/BGVAR")
# library(readxl)
# library(tidyverse)
# library(BGVAR)
# library(glue)
# options(scipen=9999)
# path_raw <- '../../Science/Regional-Multipliers-RU/Data/РегДанные.xlsx'
# path_data <- '../../Science/Regional-Multipliers-RU/Results/ModelData/'
# path_result <- '../../Science/Regional-Multipliers-RU/Results/'
# flag <- ''
# saver <- function(object, folder, name) {
#   saveRDS(object, paste0(folder, '/', name), compress = FALSE)
#   gc()
# }
# gamma_coef <- function(mode, sd) {
#   # Returns a list with shape k and scale parameter theta
#   mode_sq <- mode ^ 2
#   sd_sq <- sd ^ 2
#   k <- (2 + mode_sq / sd_sq + sqrt((4 + mode_sq / sd_sq) * mode_sq / sd_sq)) / 2
#   theta <- sqrt(sd_sq / k)
#   return(list("k" = k, "theta" = theta))
# }
# estim_model <- function(nvar, data, W, draws, burnin, plag, prior, prmean, flag, SV, thin, OE.weights, ncores, path_result, a_1 = 0.01, b_1 = 0.01) {
#   mname <- glue("{nvar}v{paste0(plag,collapse='_')}p-{prior}-{prmean}-{ifelse(SV, 'SV', 'NSV')}-{draws}-{burnin}-{thin}{flag}")
#   path_model <- glue('{path_result}{mname}/')
#   dir.create(path_model, showWarnings = FALSE)
#   if (file.exists(glue('{path_model}{mname}.RDS'))) {
#     return(NULL)
#   } else {
#     model.1 <- bgvar(Data=data,
#                      W = W,
#                      draws = draws,
#                      burnin = burnin,
#                      plag = plag,
#                      prior = prior,
#                      hyperpara = list(prmean = prmean,
#                                       a_1 = a_1,
#                                       b_1 = b_1
#                      ),
#                      SV = SV,
#                      thin = thin,
#                      trend = T,
#                      hold.out = 0,
#                      eigen = 1.05,
#                      expert = list(OE.weights=OE.weights,
#                                    cores = ncores
#                      )
#     )
#     saver(model.1, path_model, glue('{mname}.RDS'))
#     msumm <- summary(model.1)
#     saver(list(
#       CD = msumm$CD$perc,
#       res = msumm$res$p.res,
#       crosscorr = msumm$cross.corr$res.res,
#       crosscorrdata = msumm$cross.corr$dat.res
#     ), path_model, glue('{mname}-summary.RDS'))
#   }
# }
# calc_irfs <- function(data, mname, path_result, QQ, ncores) {
#   path_model <- glue('{path_result}{mname}/')
#   model.1 <- readRDS(glue('{path_model}{mname}.RDS'))
#   ### Federal Budget
#   ident <- 'chol'
#   shockinfo <- get_shockinfo(ident)
#   shockinfo$shock <- "RU.expfed"
#   shockinfo$scale <- 1
#   if (!file.exists(glue('{path_model}{mname}-{shockinfo$shock}_{ident}.RDS'))) {
#     i3 <- irf(model.1, n.ahead=24, shockinfo=shockinfo, quantiles = QQ,
#               expert = list(use_R = F, cores = ncores, calc_median = FALSE))
#     # plot(i3, cumulative = F, resp = 'y', quantiles = QQ)
#     saver(i3, path_model, glue('{mname}-{shockinfo$shock}_{ident}.RDS'))
#   }
#   ### MP
#   ident <- 'chol'
#   shockinfo <- get_shockinfo(ident)
#   shockinfo$shock <- "RU.i"
#   shockinfo$scale <- 1
#   if (!file.exists(glue('{path_model}{mname}-{shockinfo$shock}_{ident}.RDS'))) {
#     i3 <- irf(model.1, n.ahead=24, shockinfo=shockinfo, quantiles = QQ,
#               expert = list(use_R = F, cores = ncores, calc_median = FALSE))
#     # plot(i3, cumulative = F, resp = 'y', quantiles = QQ)
#     saver(i3, path_model, glue('{mname}-{shockinfo$shock}_{ident}.RDS'))
#   }
#   ### Oil Prices
#   ident <- 'chol'
#   shockinfo <- get_shockinfo(ident)
#   shockinfo$shock <- "RU.poil"
#   shockinfo$scale <- 1
#   if (!file.exists(glue('{path_model}{mname}-{shockinfo$shock}_{ident}.RDS'))) {
#     i3 <- irf(model.1, n.ahead=24, shockinfo=shockinfo, quantiles = QQ,
#               expert = list(use_R = F, cores = ncores, calc_median = FALSE))
#     # plot(i3, cumulative = F, resp = 'y', quantiles = QQ)
#     saver(i3, path_model, glue('{mname}-{shockinfo$shock}_{ident}.RDS'))
#   }
#   ### Regional Expenditure Shocks
#   shocks <- paste0(names(data) %>% head(-1), '.exp')
#   ident <- 'chol'
#   shocks <- shocks[!shocks %in% str_sub(list.files(glue('{path_model}regions_chol')), -15, -10)]
#   dir.create(glue('{path_model}regions_chol'), showWarnings = FALSE)
#   for (i in shocks) {
#     print(i)
#     shockinfo <- get_shockinfo(ident)
#     shockinfo$shock <- i
#     shockinfo$scale <- 1
#     ii <- irf(model.1, n.ahead=24, shockinfo=shockinfo, quantiles = QQ, expert = list(cores = ncores, use_R = F, calc_median = FALSE))
#     # plot(ii, cumulative = T, resp = 'MW', quantiles = QQ)
#     gc()
#     saver(ii, glue('{path_model}regions_chol'), glue('{mname}-{shockinfo$shock}_{ident}.RDS'))
#     rm(ii)
#   }
#   if (!file.exists(glue('{path_model}{mname}-regIRFs-chol.RDS'))) {
#     f <- list.files(glue('{path_model}regions_chol'), full.names = T)
#     for (i in 1:length(f)) {
#       est <- readRDS(f[i])
#       effects <- plyr::adply(est$posterior, .margins = 1:3, .id = c('response','h','shock')) %>%
#         separate(response, into = c('region', 'variable'), sep = '\\.') %>%
#         group_by(region, variable) %>%
#         arrange(h)
#       if (i == 1) result <- effects else result <- rbind(result, effects)
#     }
#     result <- result %>%
#       separate(shock, c("source", "shock"))
#     saver(result, path_model, glue('{mname}-regIRFs-chol.RDS'))
#     gc()
#   }
# }
# calc_irfs_sign <- function(data, mname, path_result, QQ, ncores) {
#   path_model <- glue('{path_result}{mname}/')
#   model.1 <- readRDS(glue('{path_model}{mname}.RDS'))
#   ### Regional Expenditure Shocks
#   shocks <- paste0(names(data) %>% head(-2), '.exp')
#   ident <- 'sign'
#   shocks <- shocks[!shocks %in% str_sub(list.files(glue('{path_model}regions_sign')), -15, -10)]
#   dir.create(glue('{path_model}regions_sign'), showWarnings = FALSE)
#   for (i in shocks) {
#     print(i)
#     reg <- substr(i, 1, 2)
#     shockinfo <- get_shockinfo(ident)
#     # AD shock
#     shockinfo <- add_shockinfo(shockinfo, shock=paste0(reg, '.y'),
#                                restriction=c(paste0(reg,'.pi'), paste0(reg, '.exp')),
#                                sign=c(">",0), horizon=1, scale=1)
#     # AS shock
#     shockinfo <- add_shockinfo(shockinfo, shock=paste0(reg,'.pi'),
#                                restriction=c(paste0(reg,'.y'), paste0(reg, '.exp')),
#                                sign=c("<",0), horizon=1, scale=1)
#     # Fiscal shock
#     shockinfo <- add_shockinfo(shockinfo, shock=i,
#                                restriction=paste0(reg, '.y'), sign=">", horizon=3, scale=1)
#     ii <- irf(model.1, n.ahead=24, shockinfo=shockinfo, quantiles = QQ, verbose = T,
#               expert = list(cores = ncores, use_R = F, MaxTries = 5000, calc_median = FALSE))
#     gc()
#     saver(ii, glue('{path_model}regions_sign'), glue('{mname}-{i}_{ident}.RDS'))
#     rm(ii)
#   }
#   if (!file.exists(glue('{path_model}{mname}-regIRFs-chol.RDS'))) {
#     f <- list.files(glue('{path_model}regions_sign'), full.names = T)
#     for (i in 1:length(f)) {
#       est <- readRDS(f[i])
#       effects <- plyr::adply(est$posterior, .margins = 1:3, .id = c('response','h','shock')) %>%
#         separate(response, into = c('region', 'variable'), sep = '\\.') %>%
#         group_by(region, variable) %>%
#         arrange(h)
#       if (i == 1) result <- effects else result <- rbind(result, effects)
#     }
#     result <- result %>%
#       separate(shock, c("source", "shock"))
#     saver(result, path_model, glue('{mname}-regIRFs-sign.RDS'))
#     gc()
#   }
# }
# calc_irfs_joint <- function(data, mname, path_result, QQ, ncores) {
#   ## Joint Fiscal Expansion ----
#   path_model <- glue('{path_result}{mname}/')
#   model.1 <- readRDS(glue('{path_model}{mname}.RDS'))
#   explanatory <- read_excel(path_raw, sheet = 'ФАКТОРЫ СВОД') %>%
#     select(region, OKATO_id, GRP) %>%
#     left_join(codes) %>%
#     filter(!is.na(region)) %>%
#     arrange(GRP)
#   poor10 <- explanatory$abbrev[1:10]
#   rich10 <- explanatory$abbrev[70:79]
#   ident <- 'chol'
#   # smallest
#   shockinfo <- data.frame(shock = paste0(poor10, '.exp'), scale = 1, global = TRUE)
#   ii <- irf(model.1, n.ahead=24, shockinfo=shockinfo, quantiles = QQ, expert = list(cores = 1, use_R = F, calc_median=FALSE))
#   effects <- plyr::adply(ii$posterior, .margins = 1:3, .id = c('response','h','shock')) %>%
#     separate(response, into = c('region', 'variable'), sep = '\\.') %>%
#     group_by(region, variable) %>%
#     arrange(h)
#   saver(effects, path, paste0(mname, '-', 'joint_10poorest_', ident, '.RDS'))
#   # largest
#   shockinfo <- data.frame(shock = paste0(rich10, '.exp'), scale = 1, global = TRUE)
#   ii <- irf(model.1, n.ahead=24, shockinfo=shockinfo, expert = list(cores = 30, use_R = T))
#   effects <- plyr::adply(ii$posterior, .margins = 1:3, .id = c('response','h','shock')) %>%
#     separate(response, into = c('region', 'variable'), sep = '\\.') %>%
#     group_by(region, variable) %>%
#     arrange(h)
#   saver(effects, path, paste0(mname, '-', 'joint_10richest_', ident, '.RDS'))
# }
# codes <- read_excel(path_raw, sheet = 'КОДЫ') %>%
#   select(region=name_official, abbrev=abbrev, OKATO_id) %>%
#   arrange(abbrev)
# # Preliminaries ----
# ### Weighing Matrix ----
# W <- readRDS(glue('{path_data}W.RDS'))
# ### GRP Weights for Aggregation ----
# aggW <- read_excel(path_raw, sheet='ФАКТОРЫ СВОД') %>%
#   left_join(codes) %>%
#   select(abbrev, w = GRPShareSample) %>%
#   arrange(abbrev)
# ### Data ----
# data <- readRDS(glue('{path_data}data.RDS'))
# ### Check Order  ----
# all(colnames(W) == names(data[-80]))
# all(rownames(W) == names(data[-80]))
# all(colnames(W) == aggW[,1])
# ## Dominant Unit Settings ----
# RU.weights <- list()
# RU.weights$weights <- aggW %>% pull(w)
# names(RU.weights$weights) <- aggW %>% pull(abbrev)
# RU.weights$variables <- c('poil', 'expfed', 'e', 'revfed', 'i', 'pi', 'y')
# RU.weights$exo <- c('poil', 'expfed', 'e', 'revfed', 'i')
# OE.weights <- list(RU=RU.weights)
# hpp <- gamma_coef(0.2, 0.4)
# # Estimation ----
# ## MN: 50K+100K, thin=10, 1 lag ----
# nvar <- 3
# draws <- 100
# burnin <- 100
# thin = 10
# plag <- c(1,1)
# prior <- 'MN'
# prmean <- 0.9
# SV <- TRUE
# ncores <- 1
# QQ <- c(0.16, 0.5, 0.84)
# mname <- glue("{nvar}v{paste0(plag,collapse='_')}p-{prior}-{prmean}-{ifelse(SV, 'SV', 'NSV')}-{draws}-{burnin}-{thin}{flag}")
# 
# path_model <- glue('{path_result}{mname}/')
# model.1 <- readRDS(glue('{path_model}{mname}.RDS'))
# explanatory <- read_excel(path_raw, sheet = 'ФАКТОРЫ СВОД') %>%
#   select(region, OKATO_id, GRP) %>%
#   left_join(codes) %>%
#   filter(!is.na(region)) %>%
#   arrange(GRP)
# poor10 <- explanatory$abbrev[1:10]
# rich10 <- explanatory$abbrev[70:79]
# ident <- 'chol'
# # smallest
# shockinfo <- data.frame(shock = paste0(poor10, '.exp'), scale = 1, global = TRUE)
# ii <- irf(model.1, n.ahead=24, shockinfo=shockinfo, quantiles = QQ, expert = list(cores = ncores, use_R = F, calc_median=FALSE))
# 
# calc_irfs_joint(