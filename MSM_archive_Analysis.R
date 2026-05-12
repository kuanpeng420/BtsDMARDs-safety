  # Marginal Structural Model: Comparative safety of b/ts DMARDs in RA
  # Data source: CDARS, Hong Kong Hospital Authority; cohort 2009-2022, follow-up to 2023-12-31
  #MSM----
  pacman::p_load(tidyverse,data.table,readxl,writexl,broom,lubridate,lmtest,sandwich,survival,mice,boot,tictoc,doParallel,tableone,parallel,nnet,ggplot2,gridExtra,grid,ggpubr)
  ori_trial <- readRDS("msm_data_20250420_corrected.RDS")
  dx <- readRDS("Diagnosis.RDS")
  dx$date <- as_date(dx$date)
  setorder(dx,date) 
  ori_trial$DateofRegisteredDeath <- as_date(ori_trial$DateofRegisteredDeath)
  outcome <- data.table(read_xlsx("Outcome_list.xlsx"))
  # Compute Charlson Comorbidity Index score from individual comorbidity indicators
  numcol <- colnames(ori_trial)[str_detect(colnames(ori_trial),"cci")]
  ori_trial[,(numcol):= lapply(.SD, as.character),.SDcol=numcol][, (numcol):=lapply(.SD, as.numeric),.SDcol=numcol]
  ori_trial[, cci := cci.mi+cci.chf+cci.pvd+cci.cbd+cci.copd+cci.dementia+cci.paralysis+(cci.dm_com0&!cci.dm_com1)+cci.dm_com1*2+cci.crf*2+(cci.liver_mild&!cci.liver_modsev)+cci.liver_modsev*3+cci.ulcers+cci.rd*6+(cci.cancer&!cci.cancer_mets)*2+cci.cancer_mets*6]
  # Exclude agents with <200 treatment episodes to avoid extreme variability in small samples
  underpower_agent <- ori_trial[,uniqueN(ReferenceKey),.(agent)][V1<200]$agent
  ori_trial <- ori_trial[!agent%in%underpower_agent]
  #Remove observation end date after 2023-09-30
  ori_trial <- ori_trial[!timepoint>as.Date("2023-09-30")]
  ori_trial$moa <- factor(ori_trial$moa,levels = c("TNFi","JAKi","IL6","Lymphocyte"))
  ori_trial$agent <- factor(ori_trial$agent,levels = c("ETANERCEPT","ADALIMUMAB","GOLIMUMAB","TOCILIZUMAB","TOFACITINIB","ABATACEPT","RITUXIMAB","BARICITINIB"))    
  moalist <- unique(ori_trial$moa)

  # Common covariate vectors ----
  # fv_num: baseline covariates only (numerator of stabilised IPCW)
  # fv_denom: baseline + time-varying covariates updated every 3 months (denominator of stabilised IPCW)
  # Time-varying confounders include ESR/CRP (disease activity markers) and concomitant medications
  # (csDMARDs, glucocorticoids, NSAIDs, opioids) that reflect inadequate disease control
  fv_num <- c("visit","cci","sex","age","base_btsdmards","base_csdmards",
              "prednisolone_b","nsaid_b","opioid_b","esr_b","crp_b",
              "HYDROXYCHLOROQUINE_b","SULPHASALAZINE_b","METHOTREXATE_b",
              "LEFLUNOMIDE_b","AZATHIOPRINE_b","CYCLOSPORIN_b")
  fv_denom <- c(fv_num, "prednisolone","nsaid","opioid","esr","crp",
                "HYDROXYCHLOROQUINE","SULPHASALAZINE","METHOTREXATE",
                "LEFLUNOMIDE","AZATHIOPRINE","CYCLOSPORIN")
  formula_num <- as.formula(paste("trt ~", paste(fv_num, collapse = "+")))
  formula_denom <- as.formula(paste("trt ~", paste(fv_denom, collapse = "+")))

  covars_plr <- "sex + age + cci + base_btsdmards + base_csdmards + visit + prednisolone_b + nsaid_b + opioid_b + esr_b + crp_b + HYDROXYCHLOROQUINE_b + AZATHIOPRINE_b + SULPHASALAZINE_b + METHOTREXATE_b + LEFLUNOMIDE_b + CYCLOSPORIN_b"
  covars_plr_novisit <- "cci + sex + age + base_btsdmards + base_csdmards + prednisolone_b + nsaid_b + opioid_b + esr_b + crp_b + HYDROXYCHLOROQUINE_b + SULPHASALAZINE_b + METHOTREXATE_b + LEFLUNOMIDE_b + AZATHIOPRINE_b + CYCLOSPORIN_b"
  plr_formula <- function(dv, grp) as.formula(paste(dv, "~", grp, "+", covars_plr))
  # plr_formula_ix: includes visit × agent interaction to allow non-proportional hazards over time
  plr_formula_ix <- function(dv, grp) as.formula(paste(dv, "~ visit + visit *", grp, "+", grp, "+", covars_plr_novisit))

  #Helper----
    fmt <- function(x,n) format(round(x,n), nsmall=n)
    as.data.table.glm <- function(m, robust=F, cluster=NULL, vcov.=NULL) {
      if(!is.null(vcov.)) m.vcov <- vcov.
      else if(!is.null(cluster)) m.vcov <- multiwayvcov::cluster.vcov(m, m$data[[cluster]])
      else if(robust) m.vcov <- sandwich::vcovHC(m, type="HC1")
      else m.vcov <- vcov(m)
      t <- coeftest(m, vcov.=m.vcov)
      data.table(fmt(cbind(Estimate=exp(t[,"Estimate"]), exp(coefci(m, vcov.=m.vcov)), p=t[,"Pr(>|z|)"]), 2), keep.rownames=T)
    }      
    
    ctb <- function(ob){
      ctb_vars <- c("sex", "age", "cci", "cci.mi", "cci.chf", "cci.pvd", "cci.cbd", "cci.copd",
                    "cci.dementia", "cci.rd", "cci.ulcers", "cci.liver_mild", "cci.dm_com0", "cci.dm_com1",
                    "cci.paralysis", "cci.crf", "cci.cancer", "cci.cancer_mets", "cci.liver_modsev",
                    "base_btsdmards", "base_csdmards", "prednisolone_b", "nsaid_b", "opioid_b",
                    "HYDROXYCHLOROQUINE_b", "SULPHASALAZINE_b", "METHOTREXATE_b", "LEFLUNOMIDE_b",
                    "AZATHIOPRINE_b", "CYCLOSPORIN_b", "esr_b", "crp_b")
      ctb_data <- ori_trial[visit == 0, c(ob, ctb_vars), with = FALSE]
      ctb_data[, sex := as.factor(sex)]

      print(CreateTableOne(vars = ctb_vars, data = ctb_data, strata = ob,
                           factorVars = c("sex", numcol, "HYDROXYCHLOROQUINE_b", "SULPHASALAZINE_b", "METHOTREXATE_b", "LEFLUNOMIDE_b", "AZATHIOPRINE_b", "CYCLOSPORIN_b"),
                           addOverall = T, test = T),
            showAllLevels = F, quote = F, noSpaces = T, test = T)
    }
    
    summarize_trial <- function(description, moa_col,event) {
      # Ensure the column exists
      if (!moa_col %in% colnames(trial)) stop(paste("Column", moa_col, "not found in the dataset."))
    
      # Create the event summary
      trialsum <- trial %>%
        filter(visit == maxVisit_cens) %>%
        count(!!sym(moa_col), !!sym(event)) %>%
        pivot_wider(names_from = !!sym(event), values_from = n, values_fill = 0) %>%
        mutate(`0` = `0` + `1`, Description = description)
    
      # Create the follow-up summary
      fltsum <- trial %>%
        filter(cens_new == 1) %>%
        mutate(flt = as.numeric(obenddate - start)) %>%
        group_by(!!sym(moa_col)) %>%
        summarise(
          `Median follow-up time (Days)` = paste0(
            median(flt), " (", quantile(flt, 0.25), ", ", quantile(flt, 0.75), ")"
          ),
          `Person years` = as.integer(sum(flt) / 365.25)
        )
    
      # Merge and calculate crude incidence rate
      msmsum <- trialsum %>%
        left_join(fltsum, by = moa_col) %>%
        mutate(
          `Crude incidence rate per 1000 person years` = round(`1` / `Person years` * 1000, 1)
        ) %>%
        select(-`Person years`)
    
      return(msmsum)
    }    
    
      # Standardised cumulative incidence via g-computation on the pooled logistic regression model
      # Used for all-cause mortality (no competing event)
      pseudofunction <- function(data,rgs,lis,grp,n){
      pseudo <- data %>% filter(visit == 0) %>% slice(rep(1:n(), each = 20*length(lis)))
      pseudo <- pseudo %>%
        mutate(
          visit = rep(0:19, times = n * length(lis)),
          {{grp}} := rep(rep(lis,each=20),n)
        )      
      pseudo <- pseudo %>% mutate(p=1 - predict(rgs, newdata = pseudo, type = "response")) %>% group_by(temname, {{grp}}) %>% mutate(s=cumprod(p))
      results <- rbind(tibble({{grp}}:=lis,visit=0,mean_incidence=0),pseudo %>% group_by(visit, {{grp}})%>%
        summarize(mean_incidence =  mean(1 - s),.groups = "keep")%>%
        ungroup()%>% mutate(visit = visit+1))        
      return(results)
      }
      
      
      # Discrete-time competing risks CIF: models event and death simultaneously
      # S(t) = cumulative product of P(no event, no death); CIF_event = sum[S(t-1) * P(event|t)]
      pseudo_cr_function <- function(data, rgsE, rgsD, lis, grp, T = 20L) {
      stopifnot(is.character(grp), length(grp) == 1L)

      # Variables needed by both models
      varsE <- all.vars(formula(rgsE))
      varsD <- all.vars(formula(rgsD))
      need  <- unique(c(varsE, varsD, "temname", "visit", grp))
      need  <- intersect(need, names(data))

      # Baseline rows only with needed columns
      data0 <- data %>%
        filter(visit == 0) %>%
        select(all_of(need))

      # Align grp factor levels and lis
      if (is.factor(data[[grp]])) {
        lvl <- levels(data[[grp]])
        data0[[grp]] <- factor(data0[[grp]], levels = lvl)
        lis_factor <- factor(lis, levels = lvl)
      } else {
        lis_factor <- lis
      }

      # Replicate each baseline row T times and assign visit 0..T-1
      n0 <- nrow(data0)
      if (n0 == 0L) stop("No baseline rows (visit == 0) found in data.")
      pseudo <- data0 %>%
        uncount(weights = T, .remove = FALSE, .id = "t_index") %>%
        mutate(visit = (t_index - 1L)) %>%
        select(-t_index)

      # Cartesian expand to all levels in lis for grp
      pseudo <- expand_grid(pseudo, setNames(list(lis_factor), grp))

      # Predict probabilities
      pseudo <- pseudo %>%
        mutate(
          pE = as.numeric(predict(rgsE, newdata = ., type = "response")),
          pD = as.numeric(predict(rgsD, newdata = ., type = "response"))
        )

      # Vectorized CIF within temname × grp
      pseudo <- pseudo %>%
        group_by(temname, .data[[grp]]) %>%
        arrange(visit, .by_group = TRUE) %>%
        mutate(
          q    = pmax(pmin(1 - pE - pD, 1), 0),
          S    = lag(cumprod(q), default = 1),
          incE = S * pE,
          cifE = cumsum(incE)
        ) %>%
        ungroup()

      # Average CIF across individuals and add visit=0 row
      res <- bind_rows(
        tibble(!!grp := lis_factor, visit = 0L, mean_CIF_event = 0),
        pseudo %>%
          group_by(visit, .data[[grp]]) %>%
          summarise(mean_CIF_event = mean(cifE), .groups = "drop") %>%
          mutate(visit = visit + 1L)  # shift so visit=1 corresponds to first interval end
      )

      res
    }

    # Identify the visit at which treatment deviation (switch or discontinuation) first occurs
    # Patients are censored at that point; lag of 90 days (183 for rituximab) is embedded in data prep
    setup_censoring <- function(trial, has_event = TRUE) {
      trial <- trial %>%
        mutate(trt_change = ifelse(trt == 0, 1, 0)) %>%
        group_by(ReferenceKey, agent) %>%
        mutate(
          cens_visit = ifelse(is.na((which(trt_change %in% 1)[1])), max(visit), (which(trt_change %in% 1)[1] - 2)),
          cens_new = ifelse(cens_visit == visit, 1, ifelse(cens_visit < visit, NA, 0)),
          maxVisit_cens = ifelse(cens_visit == -1, 0, cens_visit),
          deathOverall_new = ifelse(is.na(cens_new), NA,
                                    ifelse(is.na((which(death %in% 1)[1])), 0,
                                           ifelse((which(death %in% 1)[1] - 1) <= maxVisit_cens, 1, 0)))
        )
      if (has_event) {
        trial <- trial %>%
          mutate(
            eventOverall_new = ifelse(is.na(cens_new), NA,
                                      ifelse(is.na((which(event %in% 1)[1])), 0,
                                             ifelse((which(event %in% 1)[1] - 1) <= maxVisit_cens, 1, 0)))
          )
      }
      trial <- trial %>% ungroup()
      return(trial)
    }

    # Construct stabilised weights = IPTW (baseline, multinomial) × IPCW (time-varying, per-arm binomial)
    # Truncated at 1st/99th percentiles to reduce influence of extreme weights
    build_ipw_weights <- function(trial, agentlist) {
      trial_uncens <- trial %>% filter(visit > 0)
      nFit_list <- list()
      dFit_list <- list()
      for (i in seq_along(agentlist)) {
        agent_data <- trial_uncens[trial_uncens$agent == agentlist[i], ]
        agent_data <- droplevels(agent_data) 
        nFit_list[[tolower(agentlist[i])]] <- glm(formula_num, data = agent_data, family = binomial())
        dFit_list[[tolower(agentlist[i])]] <- glm(formula_denom, data = agent_data, family = binomial())
      }
      for (i in seq_along(agentlist)) {
        ag <- tolower(agentlist[i])
        trial[[paste0("pnum_", ag)]] <- predict(nFit_list[[ag]], newdata = trial, type = 'response')
        trial[[paste0("pdenom_", ag)]] <- predict(dFit_list[[ag]], newdata = trial, type = 'response')
      }
      trial <- trial[order(trial$ReferenceKey, trial$visit), ]
      ag_lower <- tolower(as.character(trial$agent))
      pnum_match <- numeric(nrow(trial))
      pdenom_match <- numeric(nrow(trial))
      for (ag in unique(ag_lower)) {
        idx <- ag_lower == ag
        pnum_match[idx] <- trial[[paste0("pnum_", ag)]][idx]
        pdenom_match[idx] <- trial[[paste0("pdenom_", ag)]][idx]
      }
      trial <- trial %>%
        mutate(
          numCont = ifelse(visit == 0, 1, trt * pnum_match + (1 - trt) * (1 - pnum_match)),
          denCont = ifelse(visit == 0, 1, trt * pdenom_match + (1 - trt) * (1 - pdenom_match))
        ) %>%
        group_by(ReferenceKey, agent) %>%
        mutate(
          k1_0 = cumprod(numCont),
          k1_w = cumprod(denCont)
        ) %>%
        ungroup() 

  baseline <- trial %>% filter(visit == 0) %>%
    distinct(ReferenceKey, agent, .keep_all = TRUE)

  iptw_denom_formula <- as.formula(paste("agent ~", 
    "sex + age + cci + base_btsdmards + base_csdmards +",
    "prednisolone_b + nsaid_b + opioid_b + esr_b + crp_b +",
    "HYDROXYCHLOROQUINE_b + SULPHASALAZINE_b + METHOTREXATE_b +",
    "LEFLUNOMIDE_b + AZATHIOPRINE_b + CYCLOSPORIN_b"))

  mod_iptw_den <- multinom(iptw_denom_formula, data = baseline, trace = FALSE)
  mod_iptw_num <- multinom(agent ~ 1,          data = baseline, trace = FALSE)

  p_den_mat <- predict(mod_iptw_den, newdata = baseline, type = "probs")
  p_num_mat <- predict(mod_iptw_num, newdata = baseline, type = "probs")       
        
pick_prob <- function(mat, ag_vec) {
    mat[cbind(seq_len(nrow(mat)), match(as.character(ag_vec), colnames(mat)))]
  }
  baseline$iptw_num <- pick_prob(p_num_mat, baseline$agent)
  baseline$iptw_den <- pick_prob(p_den_mat, baseline$agent)
  baseline$iptw_sw  <- baseline$iptw_num / baseline$iptw_den        
  
  trial <- trial %>%
    left_join(baseline %>% select(ReferenceKey, agent, iptw_sw),
                     by = c("ReferenceKey", "agent"))        

  trial <- trial %>%
    mutate(
      stabw_ipcw  = k1_0 / k1_w,            # stabilised IPCW: cumulative numerator / cumulative denominator
      stabw_final = iptw_sw * stabw_ipcw     # final stabilised weight = IPTW × IPCW
    ) %>%
    filter(timepoint < obenddate)        
    trial$stabw_t <- trial$stabw_final
    trial$stabw_t[trial$stabw_final > quantile(trial$stabw_final, 0.99)] <- quantile(trial$stabw_final, 0.99)
    trial$stabw_t[trial$stabw_final < quantile(trial$stabw_final, 0.01)] <- quantile(trial$stabw_final, 0.01)
    trial$stabw_t <- round(trial$stabw_t, 4)
    return(trial)
    }

    extract_hr <- function(model, grp_var, description) {
      res <- as.data.table(model, robust = T)
      res[, `:=`(`Hazard ratio` = paste0(Estimate, " (", `2.5 %`, ",", `97.5 %`, ")"), Description = description)]
      res[[grp_var]] <- str_remove_all(res$rn, grp_var)
      res <- res[str_detect(rn, grp_var)][, .SD, .SDcols = c("Description", grp_var, "Hazard ratio")]
      return(res)
    }

    run_boot_and_collect <- function(indt, boot_fn, strata_var, R, description, estnames, conf_bonf) {
      n_workers <- max(1L, detectCores() - 4L)
      cl <- makePSOCKcluster(n_workers)
      registerDoParallel(cl)
      on.exit({ stopCluster(cl); registerDoSEQ() }, add = TRUE)
      # Load packages on worker nodes
      clusterEvalQ(cl, {
        library(dplyr); library(tidyr); library(tibble); library(data.table)
        library(reshape2); library(stats)
        NULL
      })
      # Export helper functions from the defining environment of run_boot_and_collect
      helper_source_env <- environment(run_boot_and_collect)
      helper_names <- c("covars_plr_novisit", "plr_formula_ix",
                        "pseudofunction", "pseudo_cr_function")
      helper_names <- helper_names[helper_names %in% ls(envir = helper_source_env, all.names = TRUE)]
      if (length(helper_names) > 0L) {
        clusterExport(cl, varlist = helper_names, envir = helper_source_env)
      }
      # Export all variables captured in boot_fn's closure (e.g. trial_cr, agentlist, etc.)
      boot_source_env <- environment(boot_fn)
      boot_names <- ls(envir = boot_source_env, all.names = TRUE)
      if (length(boot_names) > 0L) {
        clusterExport(cl, varlist = boot_names, envir = boot_source_env)
      }
      tic()
      set.seed(123)
      clusterSetRNGStream(cl, iseed = 123)
      reps <- boot(
        data      = indt,
        statistic = boot_fn,
        R         = R,
        strata    = indt[[strata_var]],
        parallel  = "snow",
        ncpus     = n_workers,
        cl        = cl
      )
      toc()
      # Collect point estimate + 95% CI + Bonferroni-adjusted CI
      est_dt <- data.table()
      for (i in seq_along(reps$t0)) {
        ci95   <- boot.ci(reps, type = "perc", conf = 0.95,      index = i)$percent
        ci_bon <- boot.ci(reps, type = "perc", conf = conf_bonf, index = i)$percent
        est <- data.table(
          Description   = description,
          Estimate      = estnames[i],
          Original      = paste0(round(reps$t0[i], 3), " (", round(ci95[, 4], 3), ", ", round(ci95[, 5], 3), ")"),
          Original_bonf = paste0(round(reps$t0[i], 3), " (", round(ci_bon[, 4], 3), ", ", round(ci_bon[, 5], 3), ")")
        )
        est_dt <- rbind(est_dt, est)
      }
      return(est_dt)
    }

    label_estimates <- function(dt, grplist, grp_col) {
      for (i in seq_along(grplist)) {
        dt[str_detect(Estimate, as.character(grplist[i])), (grp_col) := grplist[i]]
      }
      dt[, estype := ifelse(str_detect(Estimate, "RD"), "Risk difference", "Absolute risk")]
      return(dt)
    }

    build_final_table <- function(msmtbl, restbl, esttbl, by_col) {
      esttidy <- dcast(esttbl, formula = as.formula(paste("Description +", by_col, "~ estype")),
                       value.var = c("Original", "Original_bonf"))
      merge(merge(msmtbl, restbl, by = c("Description", by_col), allow.cartesian = T, all.x = T),
            esttidy, all.x = T)
    }

    # Bonferroni correction for 9 simultaneous outcomes (including fracture as negative control)
    conf_level_bonf <- 1 - 0.05/9 
        
   # write_xlsx(list(data.table(ctb("moa"),keep.rownames=T),data.table(ctb("agent"),keep.rownames=T)),"MSM_Tableone1.xlsx")

    #Main outcomes-------
    # Loop over outcomes defined in Outcome_list.xlsx (ICD-9-CM codes)
    # Prevalent cases within 3-year lookback are excluded on an outcome-specific basis
    msmtable_agent <- data.table()
    esttable_agent <- data.table()
    restable_agent <- data.table()
    msmtable_moa <- data.table()
    esttable_moa <- data.table()
    restable_moa <- data.table()    
    for (j in 1:length(outcome$Description)) {
    trial <- copy(ori_trial)
    event <- merge(trial,dx[str_detect(code,outcome$`ICD-9 code`[j])],all.x=T,allow.cartesian = T,by=c("ReferenceKey"))
    ongo_event <- distinct(event[date <=start & date>=start %m-% months(36)][,.(ReferenceKey,start,agent)])
    event <- event[date<=outcome_end&date>timepoint,head(.SD,1),.(ReferenceKey,start,agent)]
    trial <- merge(trial,event[,.(ReferenceKey,start,agent,code,event_date=date)],all.x=T,allow.cartesian = T,by=c("ReferenceKey","start","agent"))
    trial <- trial[!ongo_event, on = .(ReferenceKey, start, agent)]
    trial[event_date<=outcome_end&event_date>timepoint,event:=1,.(ReferenceKey,start,agent)][is.na(event),event:=0]
    trial <- copy(trial)[,obenddate:=pmin(event_date,DateofRegisteredDeath,end,start %m+% months(60),na.rm = T),c("ReferenceKey","start","agent")]
    trial[,trt:=ifelse(timepoint<obenddate,as.numeric(1),as.numeric(0))]    
    trial[event==1&event_date>end&trt==1,obenddate:=event_date]#within one period, the event_date later than the prescription end date, then the true end date should be event_date
    trial <- anti_join(trial,trial[temname%in%trial[event==1]$temname][timepoint>event_date])
    #Data set up 
    trial <- setup_censoring(trial, has_event = TRUE)
    agentlist <- unique(trial$agent)
    trial <- build_ipw_weights(trial, agentlist)
    
    msmtable_moa <- rbind(summarize_trial(outcome$Description[j],moa_col = "moa",event="event"),msmtable_moa)   
    msmtable_agent <- rbind(summarize_trial(outcome$Description[j],moa_col = "agent",event="event"),msmtable_agent)  
    
    # Competing risks: death is a competing event for all clinical outcomes
    # If event and death occur in the same interval, prioritise event
    trial_cr <- trial %>%
      mutate(
        death = if_else(event == 1L & death == 1L, 0L, death)
      )    
    
    # Pooled logistic regression without visit interaction (for odds ratio summary table)
    plrFit_SWT_moa <- glm(plr_formula("event","moa"), data=trial_cr, weights = stabw_t, family=quasibinomial())
    plrFit_SWT_agent <- glm(plr_formula("event","agent"), data=trial_cr, weights = stabw_t, family=quasibinomial())
    
    restable_moa <- rbind(restable_moa, extract_hr(plrFit_SWT_moa, "moa", outcome$Description[j]))
    restable_agent <- rbind(restable_agent, extract_hr(plrFit_SWT_agent, "agent", outcome$Description[j]))    
        
      
    mod_event_moa <- glm(plr_formula_ix("event","moa"), data = trial_cr, family = quasibinomial(), weights = stabw_t)
    mod_event_agent <- glm(plr_formula_ix("event","agent"), data = trial_cr, family = quasibinomial(), weights = stabw_t)
    mod_death_moa <- glm(plr_formula_ix("death","moa"), data = trial_cr, family = quasibinomial(), weights = stabw_t)
    mod_death_agent <- glm(plr_formula_ix("death","agent"), data = trial_cr, family = quasibinomial(), weights = stabw_t)
    

   results_moa_cr <- pseudo_cr_function(data=trial_cr,rgsE=mod_event_moa,rgsD=mod_death_moa, lis=moalist,grp="moa") %>% mutate(Month=seq(0,60,3),.by=moa)
   results_agent_cr <- pseudo_cr_function(data=trial_cr,rgsE=mod_event_agent,rgsD=mod_death_agent, lis=agentlist,grp="agent") %>% mutate(Month=seq(0,60,3),.by=agent)

      assign(paste0("pm",j),ggplot(results_moa_cr, aes(x=Month, y=mean_CIF_event)) + 
      geom_line(aes(colour=moa),linetype=2,linewidth=1) +
      geom_point(aes(colour=moa),size=3)+
      scale_x_continuous(limits = c(0, 60), breaks=seq(0,60,12)) +
      theme( axis.text = element_text(size = 20)) +  
      scale_color_brewer(palette = "Dark2") +
      ggtitle(outcome$Description[j]) +
      labs(colour="B/ts DMARDs mechanism of action") + theme_classic2() +
      theme(legend.position="none",
      axis.text = element_text(size = 18),
      title=element_text(size = 18),
      axis.title = element_blank()))
      
      assign(paste0("pa",j),ggplot(results_agent_cr %>% mutate(agent = str_to_title(agent)), aes(x=Month, y=mean_CIF_event)) + 
      geom_line(aes(colour=agent),linetype=2,linewidth=1) +
      geom_point(aes(colour=agent),size=3)+
      scale_x_continuous(limits = c(0, 60), breaks=seq(0,60,12)) +
      theme( axis.text = element_text(size = 20)) +  
      scale_color_brewer(palette = "Dark2") +
      ggtitle(outcome$Description[j]) +
      labs(colour="B/ts DMARDs agents") + theme_classic2() +
      theme(legend.position="none",
      axis.text = element_text(size = 18),
      title=element_text(size = 18),
      axis.title = element_blank()))
      
    # Stratified bootstrap (1000 replicates) to obtain percentile CIs
    # Bonferroni-adjusted CI level applied in run_boot_and_collect
    indt <- trial_cr %>%filter(!is.na(eventOverall_new)) %>%
    distinct(ReferenceKey, agent, eventOverall_new)   
    
    rsq_function_moa <- function( data, indices) {
      
      index <- data[indices,]  
      index <- index %>% select(-eventOverall_new)
      trial_s <- semi_join(trial_cr,index,by=c("ReferenceKey","agent"))
      
      plrFit_boots_moa_event <- glm(plr_formula_ix("event","moa"), data=trial_s, weights = stabw_t, family=quasibinomial())
      plrFit_boots_moa_death <- glm(plr_formula_ix("death","moa"), data=trial_s, weights = stabw_t, family=quasibinomial())
      results <- pseudo_cr_function(data=trial_s,rgsE=plrFit_boots_moa_event,rgsD=plrFit_boots_moa_death, lis=moalist,grp="moa")
      wideres <- dcast(results, visit ~ moa, value.var = 'mean_CIF_event')
      #Risk difference
      wideres$RD_JAKi <- wideres$JAKi - wideres$TNFi
      wideres$RD_IL6 <- wideres$IL6 - wideres$TNFi
      wideres$RD_Lymphocyte <- wideres$Lymphocyte - wideres$TNFi

      wideres_end <- wideres %>% filter(visit==20)
      return(c(wideres_end$TNFi,wideres_end$JAKi,wideres_end$IL6,wideres_end$Lymphocyte,
               wideres_end$RD_JAKi,wideres_end$RD_IL6,wideres_end$RD_Lymphocyte))
    }      

    estname_moa <- c("TNFi","JAKi","IL6","Lymphocyte","RD_JAKi","RD_IL6","RD_Lymphocyte")
    esttable_moa <- rbind(esttable_moa, run_boot_and_collect(indt, rsq_function_moa, "eventOverall_new", R=1000, outcome$Description[j], estname_moa, conf_level_bonf))  
       
    rsq_function_agent <- function( data, indices) {
    
    index <- data[indices,]  
    index <- index %>% select(-eventOverall_new)
    trial_s <- semi_join(trial_cr,index,by=c("ReferenceKey","agent"))
    
    plrFit_boots_agent_event <- glm(plr_formula_ix("event","agent"), data=trial_s, weights = stabw_t, family=quasibinomial())
    plrFit_boots_agent_death <- glm(plr_formula_ix("death","agent"), data=trial_s, weights = stabw_t, family=quasibinomial())
    
    results <- pseudo_cr_function(data=trial_s,rgsE=plrFit_boots_agent_event,rgsD=plrFit_boots_agent_death,lis=agentlist,grp="agent")
    wideres <- dcast(results, visit ~ agent, value.var = 'mean_CIF_event')
    #Risk difference
    wideres$RD_ADALIMUMAB <- wideres$ADALIMUMAB - wideres$ETANERCEPT
    wideres$RD_GOLIMUMAB <- wideres$GOLIMUMAB - wideres$ETANERCEPT
    wideres$RD_TOCILIZUMAB <- wideres$TOCILIZUMAB - wideres$ETANERCEPT
    wideres$RD_TOFACITINIB <- wideres$TOFACITINIB - wideres$ETANERCEPT
    wideres$RD_ABATACEPT <- wideres$ABATACEPT - wideres$ETANERCEPT
    wideres$RD_RITUXIMAB <- wideres$RITUXIMAB - wideres$ETANERCEPT
    wideres$RD_BARICITINIB <- wideres$BARICITINIB - wideres$ETANERCEPT
    wideres_end <- wideres %>% filter(visit==20)
    return(c(wideres_end$ETANERCEPT,wideres_end$ADALIMUMAB,wideres_end$GOLIMUMAB,wideres_end$TOCILIZUMAB,
             wideres_end$TOFACITINIB,wideres_end$ABATACEPT,wideres_end$RITUXIMAB,wideres_end$BARICITINIB,
             wideres_end$RD_ADALIMUMAB,wideres_end$RD_GOLIMUMAB,wideres_end$RD_TOCILIZUMAB,
             wideres_end$RD_TOFACITINIB,wideres_end$RD_ABATACEPT,wideres_end$RD_RITUXIMAB,wideres_end$RD_BARICITINIB))
  }
    
    estname_agent <- c("ETANERCEPT","ADALIMUMAB","GOLIMUMAB","TOCILIZUMAB","TOFACITINIB","ABATACEPT","RITUXIMAB","BARICITINIB","RD_ADALIMUMAB","RD_GOLIMUMAB","RD_TOCILIZUMAB","RD_TOFACITINIB","RD_ABATACEPT","RD_RITUXIMAB","RD_BARICITINIB")
    esttable_agent <- rbind(esttable_agent, run_boot_and_collect(indt, rsq_function_agent, "eventOverall_new", R=1000, outcome$Description[j], estname_agent, conf_level_bonf))
    
    }
  
    #Death---- 
    # All-cause mortality: no competing event; pseudofunction (not pseudo_cr_function) is used
    trial <- copy(ori_trial)
    trial <- copy(trial)[,obenddate:=pmin(DateofRegisteredDeath,end,start %m+% months(60),na.rm = T),c("ReferenceKey","start","agent")]
    trial[,trt:=ifelse(timepoint<obenddate,as.numeric(1),as.numeric(0))]    
    trial[death==1&DateofRegisteredDeath>end&trt==1,obenddate:=DateofRegisteredDeath]
    trial <- anti_join(trial,trial[temname%in%trial[death==1]$temname][timepoint>DateofRegisteredDeath])    
    trial <- setup_censoring(trial, has_event = FALSE)
    agentlist <- unique(trial$agent)
    trial <- build_ipw_weights(trial, agentlist)

    
    
    msmtable_moa <- rbind(summarize_trial(description="Death",moa_col = "moa",event="death"),msmtable_moa)   
    msmtable_agent <- rbind(summarize_trial("Death",moa_col = "agent",event="death"),msmtable_agent)  
    
    plrFit_SWT_moa <- glm(plr_formula("death","moa"), data=trial, weights = stabw_t, family=quasibinomial())
    plrFit_SWT_agent <- glm(plr_formula("death","agent"), data=trial, weights = stabw_t, family=quasibinomial())
    
    result_plr_moa <- extract_hr(plrFit_SWT_moa, "moa", "Death")
    result_plr_agent <- extract_hr(plrFit_SWT_agent, "agent", "Death")
    restable_moa <- rbind(restable_moa, result_plr_moa)
    restable_agent <- rbind(restable_agent, result_plr_agent)
    
    plrFit_ix_SWT_moa <- glm(plr_formula_ix("death","moa"), data=trial, weights = stabw_t, family=quasibinomial())
    plrFit_ix_SWT_agent <- glm(plr_formula_ix("death","agent"), data=trial, weights = stabw_t, family=quasibinomial())
    
    results_moa <- pseudofunction(data=trial,rgs=plrFit_ix_SWT_moa,lis=moalist,grp=moa,n=nrow(distinct(trial,ReferenceKey,agent,start))) %>% mutate(Month=seq(0,60,3),.by=moa)
    results_agent <- pseudofunction(data=trial,rgs=plrFit_ix_SWT_agent,lis=agentlist,grp=agent,n=nrow(distinct(trial,ReferenceKey,agent,start))) %>% mutate(Month=seq(0,60,3),.by=agent)
     
    pm_death <- ggplot(results_moa, aes(x=Month, y=mean_incidence)) + 
    geom_line(aes(colour=moa),linetype=2,linewidth=1) +
    geom_point(aes(colour=moa),size=3)+
    scale_x_continuous(limits = c(0, 60), breaks=seq(0,60,12)) +
    theme( axis.text = element_text(size = 20)) +  
    scale_color_brewer(palette = "Dark2") +
    ggtitle("All-cause mortality") +
    labs(colour="B/ts DMARDs mechanism of action") + theme_classic2() +
    theme(legend.position="none",
    axis.text = element_text(size = 18),
    title=element_text(size = 18),
    axis.title = element_blank()
  )    
    pa_death <- ggplot(results_agent, aes(x=Month, y=mean_incidence)) + 
    geom_line(aes(colour=agent),linetype=2,linewidth=1) +
    geom_point(aes(colour=agent),size=3)+
    scale_x_continuous(limits = c(0, 60), breaks=seq(0,60,12)) +
    theme( axis.text = element_text(size = 20)) +  
    scale_color_brewer(palette = "Dark2") +
    ggtitle("All-cause mortality") +
    labs(colour="B/ts DMARDs agents") + theme_classic2() +
    theme(legend.position="none",
    axis.text = element_text(size = 18),
    title=element_text(size = 18),
    axis.title = element_blank()
  )    
       
  
    indt <- trial %>%filter(!is.na(deathOverall_new)) %>%
    distinct(ReferenceKey, agent, deathOverall_new)      
    
    rsq_function_moa <- function( data, indices) {
      
      index <- data[indices,]  
      index <- index %>% select(-deathOverall_new)
      trial_s <- semi_join(trial,index,by=c("ReferenceKey","agent"))
      
      plrFit_boots_moa <- glm(plr_formula_ix("death","moa"), data=trial_s, weights = stabw_t, family=quasibinomial())
      results <- pseudofunction(data=trial_s,rgs=plrFit_boots_moa,lis=moalist,grp=moa,n=nrow(distinct(trial_s,ReferenceKey,agent,start)))
      wideres <- dcast(results, visit ~ moa, value.var = 'mean_incidence')      
      
      #Risk difference
      wideres$RD_JAKi <- wideres$JAKi - wideres$TNFi
      wideres$RD_IL6 <- wideres$IL6 - wideres$TNFi
      wideres$RD_Lymphocyte <- wideres$Lymphocyte - wideres$TNFi

      wideres_end <- wideres %>% filter(visit==20)
      return(c(wideres_end$TNFi,wideres_end$JAKi,wideres_end$IL6,wideres_end$Lymphocyte,
               wideres_end$RD_JAKi,wideres_end$RD_IL6,wideres_end$RD_Lymphocyte))
    }

    estname_moa <- c("TNFi","JAKi","IL6","Lymphocyte","RD_JAKi","RD_IL6","RD_Lymphocyte")
    esttable_moa <- rbind(esttable_moa, run_boot_and_collect(indt, rsq_function_moa, "deathOverall_new", R=1000, "Death", estname_moa, conf_level_bonf))
    
    rsq_function_agent <- function( data, indices) {
    
    index <- data[indices,]  
    index <- index %>% select(-deathOverall_new)
    trial_s <- semi_join(trial,index,by=c("ReferenceKey","agent"))
    
    plrFit_boots_agent <- glm(plr_formula_ix("death","agent"), data=trial_s, weights = stabw_t, family=quasibinomial())
    results <- pseudofunction(data=trial_s,rgs=plrFit_boots_agent,lis=agentlist,grp=agent,n=nrow(distinct(trial_s,ReferenceKey,agent,start)))
    wideres <- dcast(results, visit ~ agent, value.var = 'mean_incidence')
    #Risk difference
    wideres$RD_ADALIMUMAB <- wideres$ADALIMUMAB - wideres$ETANERCEPT
    wideres$RD_GOLIMUMAB <- wideres$GOLIMUMAB - wideres$ETANERCEPT
    wideres$RD_TOCILIZUMAB <- wideres$TOCILIZUMAB - wideres$ETANERCEPT
    wideres$RD_TOFACITINIB <- wideres$TOFACITINIB - wideres$ETANERCEPT
    wideres$RD_ABATACEPT <- wideres$ABATACEPT - wideres$ETANERCEPT
    wideres$RD_RITUXIMAB <- wideres$RITUXIMAB - wideres$ETANERCEPT
    wideres$RD_BARICITINIB <- wideres$BARICITINIB - wideres$ETANERCEPT
    wideres_end <- wideres %>% filter(visit==20)
    return(c(wideres_end$ETANERCEPT,wideres_end$ADALIMUMAB,wideres_end$GOLIMUMAB,wideres_end$TOCILIZUMAB,
             wideres_end$TOFACITINIB,wideres_end$ABATACEPT,wideres_end$RITUXIMAB,wideres_end$BARICITINIB,
             wideres_end$RD_ADALIMUMAB,wideres_end$RD_GOLIMUMAB,wideres_end$RD_TOCILIZUMAB,
             wideres_end$RD_TOFACITINIB,wideres_end$RD_ABATACEPT,wideres_end$RD_RITUXIMAB,wideres_end$RD_BARICITINIB))
  }
    
    estname_agent <- c("ETANERCEPT","ADALIMUMAB","GOLIMUMAB","TOCILIZUMAB","TOFACITINIB","ABATACEPT","RITUXIMAB","BARICITINIB","RD_ADALIMUMAB","RD_GOLIMUMAB","RD_TOCILIZUMAB","RD_TOFACITINIB","RD_ABATACEPT","RD_RITUXIMAB","RD_BARICITINIB")
    esttable_agent <- rbind(esttable_agent, run_boot_and_collect(indt, rsq_function_agent, "deathOverall_new", R=1000, "Death", estname_agent, conf_level_bonf))
  
    #Hospitalization----
    # Any inpatient admission (type "I"); prevalent exclusion uses 1-year lookback (vs 3-year for other outcomes)
    trial <- copy(ori_trial)
    event <- merge(trial,dx[str_detect(type,"I")&!str_detect(code,"^714.0|^V")],all.x=T,allow.cartesian = T,by=c("ReferenceKey"))
    ongo_event <- distinct(event[date <=start & date>=start %m-% months(12)][,.(ReferenceKey,start,agent)])
    event <- event[date<=outcome_end&date>timepoint,head(.SD,1),.(ReferenceKey,start,agent)]
    trial <- merge(trial,event[,.(ReferenceKey,start,agent,code,event_date=date)],all.x=T,allow.cartesian = T,by=c("ReferenceKey","start","agent"))
    trial <- trial[!ongo_event, on = .(ReferenceKey, start, agent)]
    trial[event_date<=outcome_end&event_date>timepoint,event:=1,.(ReferenceKey,start,agent)][is.na(event),event:=0]
    trial <- copy(trial)[,obenddate:=pmin(event_date,DateofRegisteredDeath,end,start %m+% months(60),na.rm = T),c("ReferenceKey","start","agent")]
    trial[,trt:=ifelse(timepoint<obenddate,as.numeric(1),as.numeric(0))]    
    trial[event==1&event_date>end&trt==1,obenddate:=event_date]#within one period, the event_date later than the prescription end date, then the true end date should be event_date
    trial <- anti_join(trial,trial[temname%in%trial[event==1]$temname][timepoint>event_date])    
    #Data set up 
    trial <- setup_censoring(trial, has_event = TRUE)
    agentlist <- unique(trial$agent)
    trial <- build_ipw_weights(trial, agentlist)
    
    msmtable_moa <- rbind(summarize_trial("Any hospitalization",moa_col = "moa",event="event"),msmtable_moa)   
    msmtable_agent <- rbind(summarize_trial("Any hospitalization",moa_col = "agent",event="event"),msmtable_agent)  

 
    #competing risk    
    trial_cr <- trial %>%
      mutate(
        # If both 1 in same interval, count as event; set death to 0
        death = if_else(event == 1L & death == 1L, 0L, death)
      )    
      
   # Pooled logistic regression model without interaction 
    plrFit_SWT_moa <- glm(plr_formula("event","moa"), data=trial_cr, weights = stabw_t, family=quasibinomial())
    plrFit_SWT_agent <- glm(plr_formula("event","agent"), data=trial_cr, weights = stabw_t, family=quasibinomial())
    
    result_plr_moa <- extract_hr(plrFit_SWT_moa, "moa", "Any hospitalization")
    result_plr_agent <- extract_hr(plrFit_SWT_agent, "agent", "Any hospitalization")
    restable_moa <- rbind(restable_moa, result_plr_moa)
    restable_agent <- rbind(restable_agent, result_plr_agent)
    
    mod_event_moa <- glm(plr_formula_ix("event","moa"), data = trial_cr, family = quasibinomial(), weights = stabw_t)
    mod_event_agent <- glm(plr_formula_ix("event","agent"), data = trial_cr, family = quasibinomial(), weights = stabw_t)
    mod_death_moa <- glm(plr_formula_ix("death","moa"), data = trial_cr, family = quasibinomial(), weights = stabw_t)
    mod_death_agent <- glm(plr_formula_ix("death","agent"), data = trial_cr, family = quasibinomial(), weights = stabw_t)


    results_moa_cr <- pseudo_cr_function(data=trial_cr,rgsE=mod_event_moa,rgsD=mod_death_moa, lis=moalist,grp="moa") %>% mutate(Month=seq(0,60,3),.by=moa)
    results_agent_cr <- pseudo_cr_function(data=trial_cr,rgsE=mod_event_agent,rgsD=mod_death_agent, lis=agentlist,grp="agent") %>% mutate(Month=seq(0,60,3),.by=agent)


    pm_hospital <- ggplot(results_moa_cr %>%
    mutate(moa = if_else(moa == "IL6", "IL6i", moa)), aes(x=Month, y=mean_CIF_event)) + 
    geom_line(aes(colour=moa),linetype=2,linewidth=1) +
    geom_point(aes(colour=moa),size=3)+
    scale_x_continuous(limits = c(0, 60), breaks=seq(0,60,12)) +
    theme( axis.text = element_text(size = 20)) +
    ggtitle("Any hospitalization") +
    scale_color_brewer(palette = "Dark2") +
    labs(colour="B/ts DMARDs mechanism of action") + theme_classic2() +
    theme(legend.position="bottom",
    legend.key.size = unit(1, "cm"),        
    legend.text = element_text(size = 18),  
    legend.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    title=element_text(size = 18),
    axis.title = element_blank()
  )
          
    pa_hospital <- ggplot(results_agent_cr %>% mutate(agent = str_to_title(agent)), aes(x=Month, y=mean_CIF_event)) + 
    geom_line(aes(colour=agent),linetype=2,linewidth=1) +
    geom_point(aes(colour=agent),size=3)+
    scale_x_continuous(limits = c(0, 60), breaks=seq(0,60,12)) +
    ggtitle("Any hospitalization") +
    scale_color_brewer(palette = "Dark2") +
    labs(colour="B/ts DMARDs agents") + theme_classic2() +
    theme(legend.position="bottom",
    legend.key.size = unit(1, "cm"),        
    legend.text = element_text(size = 18),  
    legend.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    title=element_text(size = 18),
    axis.title = element_blank()
  )
    
    indt <- trial_cr %>%filter(!is.na(eventOverall_new)) %>%
    distinct(ReferenceKey, agent, eventOverall_new)   
    
    rsq_function_moa <- function( data, indices) {
      
      index <- data[indices,]  
      index <- index %>% select(-eventOverall_new)
      trial_s <- semi_join(trial_cr,index,by=c("ReferenceKey","agent"))
      
      plrFit_boots_moa_event <- glm(plr_formula_ix("event","moa"), data=trial_s, weights = stabw_t, family=quasibinomial())
      plrFit_boots_moa_death <- glm(plr_formula_ix("death","moa"), data=trial_s, weights = stabw_t, family=quasibinomial())
      
      results <- pseudo_cr_function(data=trial_s,rgsE=plrFit_boots_moa_event,rgsD=plrFit_boots_moa_death,lis=moalist,grp="moa")
      wideres <- dcast(results, visit ~ moa, value.var = 'mean_CIF_event')
      #Risk difference
      wideres$RD_JAKi <- wideres$JAKi - wideres$TNFi
      wideres$RD_IL6 <- wideres$IL6 - wideres$TNFi
      wideres$RD_Lymphocyte <- wideres$Lymphocyte - wideres$TNFi

      wideres_end <- wideres %>% filter(visit==20)
      return(c(wideres_end$TNFi,wideres_end$JAKi,wideres_end$IL6,wideres_end$Lymphocyte,
               wideres_end$RD_JAKi,wideres_end$RD_IL6,wideres_end$RD_Lymphocyte))
    }      

    estname_moa <- c("TNFi","JAKi","IL6","Lymphocyte","RD_JAKi","RD_IL6","RD_Lymphocyte")
    esttable_moa <- rbind(esttable_moa, run_boot_and_collect(indt, rsq_function_moa, "eventOverall_new", R=1000, "Any hospitalization", estname_moa, conf_level_bonf))
        
   
    rsq_function_agent <- function( data, indices) {
    
    index <- data[indices,]  
    index <- index %>% select(-eventOverall_new)
    trial_s <- semi_join(trial_cr,index,by=c("ReferenceKey","agent"))
    
    plrFit_boots_agent_event <- glm(plr_formula_ix("event","agent"), data=trial_s, weights = stabw_t, family=quasibinomial())
    plrFit_boots_agent_death <- glm(plr_formula_ix("death","agent"), data=trial_s, weights = stabw_t, family=quasibinomial())
    
    results <- pseudo_cr_function(data=trial_s,rgsE=plrFit_boots_agent_event,rgsD=plrFit_boots_agent_death,lis=agentlist,grp="agent")
    wideres <- dcast(results, visit ~ agent, value.var = 'mean_CIF_event')
    #Risk difference
    wideres$RD_ADALIMUMAB <- wideres$ADALIMUMAB - wideres$ETANERCEPT
    wideres$RD_GOLIMUMAB <- wideres$GOLIMUMAB - wideres$ETANERCEPT
    wideres$RD_TOCILIZUMAB <- wideres$TOCILIZUMAB - wideres$ETANERCEPT
    wideres$RD_TOFACITINIB <- wideres$TOFACITINIB - wideres$ETANERCEPT
    wideres$RD_ABATACEPT <- wideres$ABATACEPT - wideres$ETANERCEPT
    wideres$RD_RITUXIMAB <- wideres$RITUXIMAB - wideres$ETANERCEPT
    wideres$RD_BARICITINIB <- wideres$BARICITINIB - wideres$ETANERCEPT
    wideres_end <- wideres %>% filter(visit==20)
    return(c(wideres_end$ETANERCEPT,wideres_end$ADALIMUMAB,wideres_end$GOLIMUMAB,wideres_end$TOCILIZUMAB,
             wideres_end$TOFACITINIB,wideres_end$ABATACEPT,wideres_end$RITUXIMAB,wideres_end$BARICITINIB,
             wideres_end$RD_ADALIMUMAB,wideres_end$RD_GOLIMUMAB,wideres_end$RD_TOCILIZUMAB,
             wideres_end$RD_TOFACITINIB,wideres_end$RD_ABATACEPT,wideres_end$RD_RITUXIMAB,wideres_end$RD_BARICITINIB))
  }
    
    estname_agent <- c("ETANERCEPT","ADALIMUMAB","GOLIMUMAB","TOCILIZUMAB","TOFACITINIB","ABATACEPT","RITUXIMAB","BARICITINIB","RD_ADALIMUMAB","RD_GOLIMUMAB","RD_TOCILIZUMAB","RD_TOFACITINIB","RD_ABATACEPT","RD_RITUXIMAB","RD_BARICITINIB")
    esttable_agent <- rbind(esttable_agent, run_boot_and_collect(indt, rsq_function_agent, "eventOverall_new", R=1000, "Any hospitalization", estname_agent, conf_level_bonf))
    
    esttable_moa <- label_estimates(esttable_moa, moalist, "moa")
    esttable_agent <- label_estimates(esttable_agent, agentlist, "agent")
    
    # Merge event counts, crude rates, weighted OR, and bootstrapped risk differences into final tables
    finaltbl <- list(build_final_table(msmtable_moa, restable_moa, esttable_moa, "moa"),
                     build_final_table(msmtable_agent, restable_agent, esttable_agent, "agent"))   
    
    write_xlsx(finaltbl,"MSM_Results_cr_IPTW+IPCW.xlsx")

    #Plot function----
    get_legend<-function(myggplot){
      tmp <- ggplot_gtable(ggplot_build(myggplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)
    }
    legendm <- get_legend(pm_hospital)   
    legenda <- get_legend(pa_hospital)   
    
    figurem <- ggarrange(pm1+theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
    pm2+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    pm3+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    # pm4+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    pm5+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    pm6+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    pm7+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    pm_death+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    pm_hospital+theme(axis.title.x = element_blank(),axis.title.y = element_blank())+ theme(legend.position="none"),
                    ncol = 4, nrow = 2,legend.grob = legendm,legend="top")    
    

    figurem <- annotate_figure(figurem,
                
                left = text_grob("Culmulative incidence of events",
                                     size = 24,rot=90),
                bottom = text_grob("Months",size = 24),
                )

    tiff("MSM_safety_moa_revision_IPTW+IPCW.jpg",width = 6800,height = 3600,res = 300,compression = "jpeg")
    figurem
    dev.off()    
    
    figurea <- ggarrange(pa1+theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
    pa2+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    pa3+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    # pa4+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    pa5+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    pa6+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    pa7+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    pa_death+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    pa_hospital+theme(axis.title.x = element_blank(),axis.title.y = element_blank())+ theme(legend.position="none"),
                    ncol = 4, nrow = 2,legend.grob = legenda,legend="top")    
    

    figurea <- annotate_figure(figurea,
                
                left = text_grob("Culmulative incidence of events",
                                     size = 24,rot=90),
                bottom = text_grob("Months",size = 24),
                )

    tiff("MSM_safety_agent_revision_IPTW+IPCW.jpg",width = 6800,height = 3600,res = 300,compression = "jpeg")
    figurea
    dev.off()      
    
