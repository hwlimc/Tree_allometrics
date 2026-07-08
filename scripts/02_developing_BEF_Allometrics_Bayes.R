getwd()
setwd('/Users/hyli0001/wrd/b/Dynamic_allometrics/')
# system('ls -alt ../')
# bp<-read.table('processed_data/plot_biomass.txt',sep='\t',head=TRUE)
# baye.var<-c("stand_id","sp_code","plot","age","d","h","Bst","Bbr","Bf","Bcr","Vst","Vst.5","ba","befa.v","befa.st","beft.v","befr.st","beft.st","bwd","PFT","Genus","Family","ftp","ft1.forest_type","ft1.dominant_prop","sdi","sdi.1","sdi.2","sdi.3","sdi_max","rsd","rsd.1","rsd.2")
# write.table(bp[,baye.var],'processed_data/plot_biomass_bayes.txt',quote=FALSE,sep='\t',row.names=FALSE)
system('ls -alt bayes_outputs')
baydata<-read.table('processed_data/plot_biomass_bayes.txt',sep='\t',head=TRUE)
library(brms)
library(posterior)
library(loo)

run_id <- format(Sys.time(), "%Y%m%d_%H%M")
out_dir <- file.path("bayes_outputs","xp_gamma_workflow", run_id)
dir.create(out_dir, recursive = TRUE, showWarnings =FALSE)

model_info <- data.frame(model = c("base_k0","ftp_k0", "ftp_k1", "sp_k0", "sp_k1","ftp_sp_k0","ftp_sp_k1","ftp_sp_k2"),
    model_file = c(
    "bayes_outputs/xp_none_k0_gamma_rsd_4-4k-99-15.rds"
    "bayes_outputs/xp_ftp_k0_gamma_rsd_4-4k-99-15.rds",
    "bayes_outputs/xp_ftp_k1_gamma_rsd_4-4k-99-15.rds",
    "bayes_outputs/xp_sp_code_k0_gamma_rsd_4-4k-99-15.rds",
    "bayes_outputs/xp_sp_code_k1_gamma_rsd_4-4k-99-15.rds",
    "bayes_outputs/xp_ftp-sp_code_k0_gamma_rsd_4-4k-99-15.rds",
    "bayes_outputs/xp_ftp-sp_code_k1_gamma_rsd_4-4k-99-15.rds",
    "bayes_outputs/xp_ftp-sp_code_k2_gamma_rsd_4-4k-99-15.rds"),
    stringsAsFactors = FALSE)


model_info$model_time<-as.character(file.info(model_info$model_file)$mtime)
model_info$model_file_size<-file.info(model_info$model_file)$size
fits<-setNames(lapply(model_info$model_file, readRDS),model_info$model)

write_txt <- function(x, file) {
    write.table(x, file.path(out_dir, file), sep = "\t",
    quote = FALSE, row.names = FALSE)}

add_meta <- function(df, model, model_file, model_time)
  {data.frame(
      run_id = run_id,
      model = model,
      model_file = model_file,
      model_time = model_time,
      df,
      row.names = NULL
    )
  }

resp_label <- function(resp) if (resp == "y1m1") "befa.st" else "befr.st"
backtransform <- function(resp, x) if (resp == "y1m1") x+1 else x
run_by_model <- function(fun) {
	do.call(rbind, Map(fun, fits, model_info$model,model_info$model_file, model_info$model_time))}

run_by_model_response <- function(fun) {
	do.call(rbind, unlist(
		Map(function(fit, model, file,time) {
      list(fun(fit, model, file, time, "y1m1"),fun(fit, model, file, time, "y2"))}, 
      fits, model_info$model, model_info$model_file,model_info$model_time),
      recursive = FALSE))}

manifest <- data.frame(
    run_id = run_id,
    run_generated_at = run_generated_at,
    model_info,
    n_rows = sapply(fits, function(x) nrow(x$data)))
write_txt(manifest, "00_model_manifest.txt")

#### STEP 1: MCMC diagnostics
get_mcmc <- function(fit, model, model_file, model_time)
  {np <- nuts_params(fit)
    s <- summarise_draws(as_draws_df(fit))
    s <- s[!grepl("^Ymi_", s$variable), ]
    out <- data.frame(
      divergences = sum(np$Parameter == "divergent__" &
      np$Value == 1),
      max_treedepth = max(np$Value[np$Parameter ==
      "treedepth__"], na.rm = TRUE),
      max_rhat = max(s$rhat, na.rm = TRUE),
      min_bulk_ess = min(s$ess_bulk, na.rm = TRUE),
      min_tail_ess = min(s$ess_tail, na.rm = TRUE))

    add_meta(out, model, model_file, model_time)}

step1 <- run_by_model(get_mcmc)
write_txt(step1, "01_mcmc_diagnostics.txt")

## STEP 2: Posterior parameter summaries

get_posteriors <- function(fit, model, model_file, model_time) {
    s <- summarise_draws(as_draws_df(fit))
    s <- s[grepl("^(b_|sd_|shape_)", s$variable), ]
    add_meta(s, model, model_file, model_time)}

step2 <- run_by_model(get_posteriors)
write_txt(step2, "02_posterior_parameter_summary.txt")

## STEP 3: Posterior predictive checks
get_ppc <- function(fit, model, model_file, model_time, resp, ndraws = 1000) {
    dat <- fit$data[!is.na(fit$data[[resp]]), , drop=FALSE]
    obs <- backtransform(resp, dat[[resp]])
    yrep <- backtransform(resp, posterior_predict(fit, newdata = dat, resp = resp, ndraws = ndraws))

    stat_names <- c("mean", "sd", "median", "q05", "q95", "max")
    stat_value <- function(x, stat) {
      switch(stat,
        mean = mean(x, na.rm = TRUE),
        sd = sd(x, na.rm = TRUE),
        median = median(x, na.rm = TRUE),
        q05 = quantile(x, 0.05, na.rm = TRUE),
        q95 = quantile(x, 0.95, na.rm = TRUE),
        max = max(x, na.rm = TRUE))}

    out <- do.call(rbind, lapply(stat_names,
    function(stat) {
      pred_stat <- apply(yrep, 1, stat_value, stat = stat)
      data.frame(
        response = resp_label(resp),
        stat = stat,
        observed = as.numeric(stat_value(obs, stat)),
        pred_median = median(pred_stat),
        pred_q05 = quantile(pred_stat, 0.05),
        pred_q95 = quantile(pred_stat, 0.95)
      )
    }))

    add_meta(out, model, model_file, model_time)
  }

step3 <- run_by_model_response(get_ppc)
write_txt(step3, "03_posterior_predictive_stats.txt")

## STEP 4: Observed versus predicted
get_obs_pred <- function(fit, model, model_file,model_time, resp) {
    dat <- fit$data[!is.na(fit$data[[resp]]), , drop=FALSE]
    pred <- fitted(fit, newdata = dat, resp = resp)[,"Estimate"]
    out <- data.frame(
      response = resp_label(resp),
      row_id = seq_len(nrow(dat)),
      x = dat$x,
      observed = backtransform(resp, dat[[resp]]),
      predicted = backtransform(resp, pred))
    out$residual <- out$observed - out$predicted
    add_meta(out, model, model_file, model_time)}

step4 <- run_by_model_response(get_obs_pred)
write_txt(step4, "04_observed_vs_predicted.txt")

pdf(file.path(out_dir, "04_observed_vs_predicted.pdf"),
width = 12, height = 7)
old_par <- par(no.readonly = TRUE)

for (resp in c("befa.st", "befr.st")) {
    par(mfrow = c(2, 4), mar = c(4, 4, 3, 1), oma = c(0,
    0, 2, 0))
    for (mod in names(fits)) {
      d <- step4[step4$response == resp & step4$model ==
      mod, ]
      lim <- range(c(d$observed, d$predicted), na.rm =
      TRUE)

      plot(d$observed, d$predicted,
        xlim = lim, ylim = lim,
        pch = 16, col = rgb(0, 0, 0, 0.35),
        xlab = "Observed", ylab = "Predicted",
        main = paste(resp, mod, sep = " | ")
      )
      abline(0, 1, col = "red", lwd = 2)
    }

    plot.new()
    mtext(paste("Observed vs predicted:", resp), outer =
    TRUE, cex = 1.2)
  }

  par(old_par)
  dev.off()

### STEP 5: PSIS-LOO and response-rsd curves
get_loo <- function(fit, model, model_file, model_time, resp) {
    dat <- fit$data[!is.na(fit$data[[resp]]), , drop=FALSE]
    lo <- loo(fit, resp = resp, newdata = dat)
    est <- lo$estimates
    pk <- pareto_k_values(lo)
    out <- data.frame(
      response = resp_label(resp),
      n = nrow(dat),
      elpd_loo = est["elpd_loo", "Estimate"],
      elpd_loo_se = est["elpd_loo", "SE"],
      p_loo = est["p_loo", "Estimate"],
      looic = est["looic", "Estimate"],
      looic_se = est["looic", "SE"],
      max_pareto_k = max(pk),
      n_pareto_k_gt_0.7 = sum(pk > 0.7),
      n_pareto_k_gt_1.0 = sum(pk > 1.0))
  add_meta(out, model, model_file, model_time)
  }

step5 <- run_by_model_response(get_loo)
step5 <- step5[order(step5$response, -step5$elpd_loo), ]
write_txt(step5, "05_loo_metrics_by_response.txt")

step5_total <- aggregate(
    cbind(elpd_loo, p_loo, looic, n, n_pareto_k_gt_0.7,n_pareto_k_gt_1.0) ~ model,
    data = step5,
    FUN = sum)
step5_total$max_pareto_k <- tapply(step5$max_pareto_k,step5$model, max)[step5_total$model]
step5_total$elpd_diff_from_best <- step5_total$elpd_loo-max(step5_total$elpd_loo)
step5_total$looic_diff_from_best <- step5_total$looic-min(step5_total$looic)
step5_total <- merge(model_info[, c("model","model_file", "model_time")], step5_total, by = "model")
step5_total<-step5_total[order(-step5_total$elpd_loo), ]
step5_total<-data.frame(run_id, run_generated_at,step5_total)

write_txt(step5_total, "05_loo_total_ranking.txt")

get_response_curve <- function(fit, model, model_file,model_time, resp, n_grid = 100) {
    dat <- fit$data[!is.na(fit$data[[resp]]), , drop =FALSE]
    x_grid <- seq(min(dat$x, na.rm = TRUE), max(dat$x,na.rm = TRUE), length.out = n_grid)
    newdat <- data.frame(x = x_grid)

    if ("h1" %in% names(fit$data)) newdat$h1 <-
    fit$data$h1[1]
    if ("h2" %in% names(fit$data)) newdat$h2 <-
    fit$data$h2[1]

    pred <- fitted(
      fit,
      newdata = newdat,
      resp = resp,
      re_formula = NA,
      probs = c(0.05, 0.95))

    obs <- data.frame(
      response = resp_label(resp),
      row_id = seq_len(nrow(dat)),
      x = dat$x,
      observed = backtransform(resp, dat[[resp]]))

    curve <- data.frame(
      response = resp_label(resp),
      x = x_grid,
      predicted = backtransform(resp, pred[, "Estimate"]),
      pred_q05 = backtransform(resp, pred[, "Q5"]),
      pred_q95 = backtransform(resp, pred[, "Q95"]))

    list(
      obs = add_meta(obs, model, model_file, model_time),
      curve = add_meta(curve, model, model_file,
      model_time))
  }

  curve_list <- unlist(Map(function(fit, model, file,
  time) {
    list(
      get_response_curve(fit, model, file, time, "y1m1"),
      get_response_curve(fit, model, file, time, "y2")
    )
  }, fits, model_info$model, model_info$model_file,
  model_info$model_time), recursive = FALSE)

  step5_obs <- do.call(rbind, lapply(curve_list, `[[`,
  "obs"))
  step5_curve <- do.call(rbind, lapply(curve_list, `[[`,
  "curve"))

write_txt(step5_obs,"05_response_vs_rsd_observed.txt")
write_txt(step5_curve,"05_response_vs_rsd_predicted.txt")

pdf(file.path(out_dir, "05_response_vs_rsd.pdf"), width=12,height=7)
  old_par <- par(no.readonly = TRUE)

  for (resp in c("befa.st", "befr.st")) {
    par(mfrow = c(2, 4), mar = c(4, 4, 3, 1), oma = c(0,
    0, 2, 0))

    for (mod in names(fits)) {
      obs_d <- step5_obs[step5_obs$response == resp &
      step5_obs$model == mod, ]
      cur_d <- step5_curve[step5_curve$response == resp &
      step5_curve$model == mod, ]
      ord <- order(cur_d$x)

      ylim <- range(c(obs_d$observed, cur_d$pred_q05,
      cur_d$pred_q95), na.rm = TRUE)

      plot(obs_d$x, obs_d$observed,
        pch = 16, col = rgb(0, 0, 0, 0.35),
        xlab = "Relative stand density (rsd)",
        ylab = resp,
        main = paste(resp, mod, sep = " | "),
        ylim = ylim
      )

      lines(cur_d$x[ord], cur_d$predicted[ord], col =
      "red", lwd = 2)
      lines(cur_d$x[ord], cur_d$pred_q05[ord], col =
      "red", lty = 2)
      lines(cur_d$x[ord], cur_d$pred_q95[ord], col =
      "red", lty = 2)
    }

    plot.new()
    mtext(paste("Response vs rsd:", resp), outer = TRUE,
    cex = 1.2)
  }

  par(old_par)
  dev.off()

  writeLines(
    c(
      paste("Run ID:", run_id),
      paste("Run generated at:", run_generated_at),
      paste("Output directory:", out_dir),
      "",
      "Files:",
      list.files(out_dir)
    ),
    file.path(out_dir, "RUN_SUMMARY.txt")
  )

  out_dir


sp_c<-readRDS("bayes_outputs/exp_decay_sp_code_rsd_4chn_4000itr_4cor_0.99del_15depth.rds")
pft_sp<-readRDS("bayes_outputs/exp_decay_PFT_sp_code_rsd_4chn_4000itr_4cor_0.99del_15depth.rds")
pft<-readRDS("bayes_outputs/exp_decay_PFT_rsd_4chn_4000itr_4cor_0.99del_15depth.rds")


ftp.sp0n<-readRDS("bayes_outputs/exp_decay_log_partial_pool_ftp_sp_code_kdepth0_rsd_student_4chn_4000itr_4cor_0.99del_15depth.rds")
ftp.sp0n
source("scripts/00_bayesian_functions.R")
pp_check(ftp.sp0n,resp='y1')
pp_check(ftp.sp0,resp='z1')


ftp.sp0
summary(ftp_sp)
summary(pft_sp)
summary(ftp)
summary(pft)
summary(sp_c)

ranef(ft_sp)$h2[, , "y1_L_Intercept"]
coef(ft_sp)



h1_cf <- coef(ft_sp)$h1
h1_param <- data.frame(
  ft1_type = dimnames(h1_cf)[[1]],
  y1_L = h1_cf[, "Estimate", "y1_L_Intercept"],
  y1_A = h1_cf[, "Estimate", "y1_A_Intercept"],
  y1_k = h1_cf[, "Estimate", "y1_k_Intercept"],
  y2_L = h1_cf[, "Estimate", "y2_L_Intercept"],
  y2_A = h1_cf[, "Estimate", "y2_A_Intercept"],
  y2_k = h1_cf[, "Estimate", "y2_k_Intercept"]
)
h1_param$y_L <- h1_param$y1_L + h1_param$y2_L
h1_param$y_x0 <- (h1_param$y1_L + h1_param$y1_A) + (h1_param$y2_L + h1_param$y2_A)
h1_param$y_A <- h1_param$y1_A + h1_param$y2_A


h2_cf <- coef(ft_sp)$h2
h2_param <- data.frame(
  sp_code = dimnames(h2_cf)[[1]],
  y1_L = h2_cf[, "Estimate", "y1_L_Intercept"],
  y1_A = h2_cf[, "Estimate", "y1_A_Intercept"],
  y1_k = h2_cf[, "Estimate", "y1_k_Intercept"],
  y2_L = h2_cf[, "Estimate", "y2_L_Intercept"],
  y2_A = h2_cf[, "Estimate", "y2_A_Intercept"],
  y2_k = h2_cf[, "Estimate", "y2_k_Intercept"]
)
h2_param$y_L <- h2_param$y1_L+h2_param$y2_L
h2_param$y_x0 <- (h2_param$y1_L+h2_param$y1_A)+(h2_param$y2_L+h2_param$y2_A)
h2_param$y_A <- h2_param$y1_A + h2_param$y2_A
head(baydata)
range(baydata$d[baydata$PFT=='ENF'])
baydata[baydata$d<0.05,]



plot(befa.st~sdi,baydata[baydata$sp_code=='PDe',],ylim=c(1,6),xlim=c(0,500))
points(befa.st~sdi,baydata[baydata$sp_code=='PD',],col=4)
points(23,5.7,col=2)
points(23,4.69,col=2)


par(mfrow=c(1,3))
x<-'rsd'
y1<-'befa.st'
y2<-'beft.st'
for (i in unique(baydata$ft1.forest_type)){
	dft <- baydata[baydata$ft1.forest_type == i, ]

plot(dft[[x]],dft[[y2]], col = 0, xlab = "Relative stand density", ylab = "Biomass expansion factors",main = i,ylim = c(0.8,4),xlim=c(0,1.1))

	sp_list<-unique(dft$sp_code)
	for (j in seq_along(sp_list)){
		sp <- sp_list[j]
		dsub<-dft[dft$sp_code==sp,]
		points(dsub[[x]], dsub[[y1]], pch = 21, bg = j,cex=0.75, col=0)
		points(dsub[[x]], dsub[[y2]], pch = 2, col = j,lwd=0.5,cex=0.75)
		}

psub <- h1_param[h1_param$ft1_type == i, ]
xseq <- seq(min(dft[[x]], na.rm = TRUE),max(dft[[x]], na.rm = TRUE),length.out = 200)
ypre1 <- psub$y1_L + psub$y1_A * exp(-psub$y1_k * xseq)
ypre2 <- psub$y_L + psub$y_A * exp(-psub$y_x0 * xseq)

lines(xseq, ypre1, col = 1, lwd = 2)
lines(xseq, ypre2, col = 2, lwd = 2, lty = 2)
}

###################### PARAMETER DISTRIBUTION SHOULD BE MADE


############################## SPECIES TEST NEEDS TO BE DONE

for (i in seq_along(sp_list)) {
  sp <- sp_list[i]

  dsub <- dft[dft$sp_code == sp, ]
  psub <- h2_param[h2_param$sp_code == sp, ]

  if (nrow(psub) == 0) next

  points(dsub[[x]], dsub[[y1]], bg = i, pch = 21)
  points(dsub[[x]], dsub[[y2]], col = i, pch = 2)

  xseq <- seq(min(dsub[[x]], na.rm = TRUE),
              max(dsub[[x]], na.rm = TRUE),
              length.out = 200)

  ypre1 <- psub$y1_L + psub$y1_A * exp(-psub$y1_k * xseq)

  ypre2 <- psub$y2_L + psub$y2_A * exp(-psub$y2_k * xseq)

  lines(xseq, ypre1, col = i, lwd = 1.5)
  lines(xseq, ypre2, col = i, lwd = 1.5,lty=2)

}

legend("topright", legend = sp_list,
       col = seq_along(sp_list), pch = 21, lwd = 1.5)


plot(befa.st~age, bppd)
points(beft.st~age, bppd[bppd$sp_code=='PDe',],col=2)
non.mod.exp<-nls(befa.st~exp_decay(rsd,L,A,k),bppd,start=c(L=11,A=12,k=3))

curve(exp_decay(x,summary(non.mod.exp)$coe[1],summary(non.mod.exp)$coe[2],summary(non.mod.exp)$coe[3]),add=TRUE)


bppd$befa.st.pre<-predict(non.mod.exp,bppd)
bppd$befa.st.res<-bppd$befa.st-bppd$befa.st.pre
curve(exp_decay(x,summary(non.mod.exp)$coe[1],summary(non.mod.exp)$coe[2],summary(non.mod.exp)$coe[3]),add=TRUE)
lines(befa.st.pre~rsd,rsd,lwd=1.1,col=4)
summary(lm(befa.st~befa.st.pre*sp_code,bppd))
plot(befa.st.res~befa.st.pre,bppd)
points(befa.st.res~befa.st.pre,bppd[bppd$sp_code=='PDe',],bg=2,pch=21)


unique(bp[,c('species','sp_code')])
head(bp)
##The number of spcies
sum(!is.na(unique(bp$sp_code)))
unique(bp$sp_code)
unique(bp$ft1.forest_type)
sum(!is.na(unique(bp$PFT)))



par(mfrow=c(2,2))
for (i in unique(bp$ft1.forest_type)){
	df1<-bp[bp$ft1.forest_type==i,]

lmx<-max(c(df1$beft.st,df1$befa.st),na.rm=TRUE)
	plot(befa.st~ rsd,df1,ylim=c(1,lmx),main=i)
	points(beft.st~rsd,df1,col=2)

	plot(befa.st~ wb.scale,ylim=c(1,lmx),df1)
	points(beft.st~ wb.scale,df1,col=2)

	plot(befa.st~sdi,ylim=c(1,lmx),df1)
	points(beft.st~ sdi,df1,col=2)

	plot(befa.st~rsd,ylim=c(1,lmx),df1)
	points(beft.st~ rsd,df1,col=2)	
	}
	

par(mfrow=c(2,2))
for (i in unique(bp$PFT)){
	df1<-bp[bp$PFT==i,]
lmx<-max(c(df1$beft.st,df1$befa.st),na.rm=TRUE)
	plot(befa.st~ wb.shape,df1,ylim=c(1,lmx),main=i)
	points(beft.st~ wb.shape,df1,col=2)

	plot(befa.st~ wb.scale,ylim=c(1,lmx),df1)
	points(beft.st~ wb.scale,df1,col=2)

	plot(befa.st~sdi,ylim=c(1,lmx),df1)
	points(beft.st~ sdi,df1,col=2)
	
	plot(befa.st~rsd,ylim=c(1,lmx),df1)
	points(beft.st~ rsd,df1,col=2)		
	}

par(mfrow=c(2,2))
for (i in unique(bp$sp_code)){
	df1<-bp[bp$sp_code==i,]
lmx<-max(c(df1$beft.st,df1$befa.st),na.rm=TRUE)
	plot(befa.st~ wb.shape,df1,ylim=c(1,lmx),main=i)
	points(beft.st~ wb.shape,df1,col=2)

	plot(befa.st~ wb.scale,ylim=c(1,lmx),df1)
	points(beft.st~ wb.scale,df1,col=2)

	plot(befa.st~sdi,ylim=c(1,lmx),df1)
	points(beft.st~ sdi,df1,col=2)
	
	plot(befa.st~rsd,ylim=c(1,lmx),df1)
	points(beft.st~ rsd,df1,col=2)		
	}



non.mod.exp<-nls(befa.st~exp_decay(rsd,L,A,k),bp[bp$PFT=='ENF',],start=c(L=1.3,A=1.4,k=.2))

plot(befa.st~rsd,bp[bp$PFT=='ENF',],col=0)
for(i in 1:7){
	sp_c<-unique(bp[bp$PFT=='ENF','sp_code'])[i]
points(befa.st~rsd,bp[bp$sp_code==sp_c&bp$PFT=='ENF',],bg=i,pch=21)}
df00<-bp[bp$PFT=='ENF',]
rsd<-as.data.frame(seq(range(df00$rsd)[1],range(df00$rsd)[2],0.001))
colnames(rsd)<-'rsd'
bp$befa.a.pre<-predict(non.mod.exp,bp)
rsd$befa.a.pre<-predict(non.mod.exp,rsd)
summary(lm(befa.a.pre~befa.st,bp))
lines(befa.a.pre~rsd,rsd,lwd=1.1,col=4)

non.mod.exp<-nls(befa.st~exp_decay(sdi,L,A,k),bp[bp$PFT=='ENF',],start=c(L=1.3,A=1.4,k=0.1))

plot(befa.st~sdi,bp[bp$PFT=='ENF',],col=0)
for(i in 1:7){
	sp_c<-unique(bp[bp$PFT=='ENF','sp_code'])[i]
points(befa.st~sdi,bp[bp$sp_code==sp_c&bp$PFT=='ENF',],bg=i,pch=21)}
df00<-bp[bp$PFT=='ENF',]
sdi<-as.data.frame(seq(range(df00$sdi)[1],range(df00$sdi)[2],0.001))
colnames(sdi)<-'sdi'
bp$befa.a.pre<-predict(non.mod.exp,bp)
sdi$befa.a.pre<-predict(non.mod.exp,sdi)
summary(lm(befa.a.pre~befa.st,bp))
lines(befa.a.pre~sdi,sdi,lwd=1.1,col=4)



unique(bp[bp$PFT=='ENF','sp_code'])[1]
plot(befa.st~sdi,bp[bp$sp_code=='PK',])
points(befa.st~sdi,bp[bp$sp_code=='PK'&bp$stand_id==1005,],col=2)


plot(befa.st~sdi,bp[bp$Family=='Pinaceae',],col=0)
for(i in 1:6){
	sp_c<-unique(bp[bp$Family=='Pinaceae','sp_code'])[i]
points(befa.st~sdi,bp[bp$sp_code==sp_c&bp$Family=='Pinaceae',],bg=i,pch=21)}
points(befa.st~sdi,bp[bp$sp_code=='LL'&bp$Family=='Pinaceae',],bg='white',pch=22,col=4)
df00<-bp[bp$Family=='Pinaceae',]
sdi<-as.data.frame(seq(range(df00$sdi)[1],range(df00$sdi)[2],1))
colnames(sdi)<-'sdi'
bp$befa.a.pre<-predict(non.mod.mm,bp)
summary(lm(befa.a.pre~befa.st,bp))
lines(befa.a.pre~sdi,sdi,lwd=1.1,col=4)


par(mfrow=c(2,3))
for (i in unique(bp$Family)){
	df1<-bp[bp$Family==i,]
lmx<-max(c(df1$beft.st,df1$befa.st),na.rm=TRUE)
	plot(befa.st~ wb.scale,ylim=c(1,lmx),df1,main=i,xlim=c(0,.6))
	points(beft.st~ wb.scale,df1,col=2)

	plot(befa.st~ wb.shape,ylim=c(1,lmx),df1,main=i,xlim=c(0,12))
	points(beft.st~ wb.shape,df1,col=2)

	plot(befa.st~sdi,ylim=c(1,lmx),df1,main=i,xlim=c(0,700))
	points(beft.st~sdi,df1,col=2)	
	}


par(mfrow=c(2,3))
for (i in unique(bp$Genus)){
	df1<-bp[bp$Genus==i,]
lmx<-max(c(df1$beft.st,df1$befa.st),na.rm=TRUE)
	plot(befa.st~ wb.scale,ylim=c(1,lmx),df1,main=i,xlim=c(0,.6))
	points(beft.st~ wb.scale,df1,col=2)

	plot(befa.st~ wb.shape,ylim=c(1,lmx),df1,main=i,xlim=c(0,12))
	points(beft.st~ wb.shape,df1,col=2)

	plot(befa.st~sdi,ylim=c(1,lmx),df1,main=i,xlim=c(0,700))
	points(beft.st~ sdi,df1,col=2)	
	
	}

par(mfrow=c(2,3))
for (i in unique(bp$sp_code)[1:22]){
	df1<-bp[bp$sp_code==i,]
lmx<-max(c(df1$beft.st,df1$befa.st),na.rm=TRUE)

	plot(befa.st~ wb.scale,ylim=c(1,lmx),df1,main=i,xlim=c(0,.6))
	points(beft.st~ wb.scale,df1,col=2)

	plot(befa.st~ wb.shape,ylim=c(1,lmx),df1,main=i,xlim=c(0,12))
	points(beft.st~ wb.shape,df1,col=2)

	plot(befa.st~sdi,ylim=c(1,lmx),df1,main=i,xlim=c(0,700))
	points(beft.st~ sdi,df1,col=2)	
	}

plot(wb.scale~sdi,bp)











library(brms)
library(cmdstanr)

source('scripts/00_bayesian_functions.R')
