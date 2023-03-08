#==============================================================================
#' get_scaling_params_Intensity_aggLevel
#' 
#' Get scaling parameters of a MRC model, where the cascade weights are 
#' considered to be dependent on the temporal aggregation level 
#' and on the intensity of precipitation that has to be disaggregated.
#' In the first step observed weights are calculated for each temporal aggregation level. 
#' Also at this stage, for each observed precipitation intensity, its index z is calculated.  
#' Second, MRC parameters are estimated by considering a dependency to 
#' temporal aggregation level, precipitation intensity class and asymmetery index class z. 
#' In the last step, scaling model are fitted on estimated MRC parameter by 
#' the method of linear least square. Plots showing estimted params and their scaling models 
#' are also returned   
#' 
#' @param vecPrecip a vector of observed precipitations 
#' @param vecDates a vector of dates, same length as \code{length(vecPrecip)}
#' @param resVecPrecip time resolution of \code{vecPrecip} in minutes 
#' @param aggLevels a vector of time aggregation levels in the cascade estimation procedure, in minutes 
#' @param by_season logical, should the parameter estimation must be done on seasonal basis. Default \code{TRUE} 
#' @param threshold for weight discarding in all estimation procedure, in mm 
#' @param threshold_alpha_int for discarding alpha estimated at small intensities, in mm
#' @param threshold_asymm for discarding weights of small intensities for estimation of asymmetry parameters, in mm
#' 
#' @return to_do a list of parameters, plots for fitting, and df of empirical MRC params
#' 
#'  
#' @author Kaltrina Maloku
#' 
#' @export
get_scaling_params_Intensity_aggLevel = function(vecPrecip,
                                                 vecDates,
                                                 resVecPrecip = 10,
                                                 by_season = T,
                                                 threshold = 0,
                                                 threshold_alpha_int = 0.8,
                                                 threshold_asymm = 0.8,
                                                 aggLevels = c(80,160,320,640,1280)){
  # estimate parameters of MRC 
  params_MRC = estimate_params_MRC(vecPrecip = vecPrecip, vecDates = vecDates,
                                   resVecPrecip = resVecPrecip,
                                   by_season = T,
                                   threshold = threshold, 
                                   threshold_alpha_int = threshold_alpha_int,
                                   threshold_asymm = threshold_asymm, 
                                   aggLevels = aggLevels)
  
  mu_start = -0.5; sigma_start = 1
  a_start = -1;b_start = -1
  nu_start = 0.5; lambda_start = 0.5
  K_start = 0.1
  
  vec.season.char = c("DJF", "MAM", "JJA", "SON")
  
  # where to predict values 
  vecIntensity = exp(seq(from = log(0.003), to = log(150), length.out = 100))
  
  # ======================= Px model =======================
  
  # First stage: estimate params of mu and sigma for each temporal scales 
  fit_intermittency_mu_sigma = params_MRC$px_aggLevel_Intensity %>% group_by(Season,aggLevel) %>% 
    dplyr::do(fit = robustbase::nlrob(value~get_Px_Intensity(Intensity = Intensity, 
                                                 mu = mu, Sigma = Sigma),
                                      control = list(maxiter = 1000, printEval = F, minFactor = 1/100000),
                          start = list(mu = mu_start, Sigma = sigma_start), maxit = 2000, data = .data)) %>%
   mutate(mu = fit$coefficients[1], Sigma = fit$coefficients[2]) %>% dplyr::select(-fit)
  
  # Second stage: estimate scaling params of mu and sigma over temporal scales, fit linear models
  fit_intermittency_mu = fit_intermittency_mu_sigma %>% group_by(Season) %>% dplyr::select(-Sigma) %>%
    dplyr::do(fit = robustbase::nlrob(mu~(a_mu*log(aggLevel/1280)+b_mu), 
                          start = list(a_mu = a_start, b_mu = b_start), 
                          data = .data, maxit = 2000)) %>%  
   mutate(a_mu = fit$coefficients[1],b_mu = fit$coefficients[2]) %>% dplyr::select(-fit)
  
  fit_intermittency_sigma = fit_intermittency_mu_sigma %>% group_by(Season) %>% dplyr::select(-mu) %>%
    dplyr::do(fit = robustbase::nlrob(Sigma~(a_sigma*log(aggLevel/1280)+b_sigma), 
                          start = list(a_sigma = a_start, b_sigma = b_start), 
                          data = .data, maxit = 2000)) %>%
   mutate(a_sigma = fit$coefficients[1],b_sigma = fit$coefficients[2]) %>%  dplyr::select(-fit)
  
  # prepare results 
  params_px = dplyr::left_join(fit_intermittency_mu,fit_intermittency_sigma, by = "Season")
  
  
  
  # ======================= Alpha model ======================= 
  # fit scaling model over temporal aggregation level 
  fit_alpha_aggLevel = params_MRC$alpha_aggLevel %>% group_by(Season) %>% 
    dplyr::do(fit = robustbase::nlrob(value~get_alpha_aggLevel(alpha0 = alpha0, H = H, 
                                                   res_aggLevel = aggLevel/1280),
                          start = list(alpha0 = 1, H = -1), maxit = 2000, data = .data)) %>%
   mutate(alpha0 = fit$coefficients[1], H = fit$coefficients[2]) %>% dplyr::select(-fit)

  # fit quadratic model in observed alpha star
  fit_quad_model = params_MRC$alpha_star_aggLevel_Intensity %>% group_by(Season) %>% 
    dplyr::do(fit = robustbase::nlrob(log(value)~(c0+c1*log(Intensity)+c2*log(Intensity)*log(Intensity)),
                          start = list(c0 = 0.1, c1 = 2, c2 = 1), 
                          data = .data, maxit = 1000,
                          control = list(maxiter = 1000, printEval = F, minFactor = 1/100000))) %>% 
   mutate(c0 = fit$coefficients[1],c1 = fit$coefficients[2], c2 = fit$coefficients[3]) %>% 
    dplyr::select(-fit) 
  
  params_alpha = dplyr::left_join(fit_alpha_aggLevel,fit_quad_model, by = "Season")
  
  
  
  # ======================= Asymmetry model =======================
  # fit linear model on estimated means for each z index class
  fit_mean_weights = params_MRC$emp_mean_z %>% dplyr::filter(!is.na(zIndex)) %>% group_by(Season) %>% 
    dplyr::do(fit = robustbase::nlrob(value~get_mean_z(z = zIndex, lambda = lambda), 
                                      control = list(maxiter = 1000, printEval = F, minFactor = 1/100000),
                          start = list(lambda = lambda_start), maxit = 2000, data = .data)) %>%
   mutate(lambda = fit$coefficients[1]) %>% dplyr::select(-fit)
  
  # fit erf model on estimated ratios for each z index class
  fit_intermittency_ratio_phi = params_MRC$phi_z %>% dplyr::filter(!is.na(zIndex)) %>% group_by(Season) %>% 
    dplyr::do(fit = robustbase::nlrob(value~get_phi_z(z = zIndex, nu = nu),
                                      control = list(maxiter = 1000, printEval = F, minFactor = 1/100000),
                          start = list(nu = nu_start), maxit = 2000, data = .data)) %>%
   mutate(nu = fit$coefficients[1]) %>% dplyr::select(-fit)
  
  params_asymm = dplyr::left_join(fit_mean_weights,fit_intermittency_ratio_phi, by = "Season")
  
  # ======================= Plot: Px model =======================
  list_plots = list()
  # predict values based on the parameters of the model
  predict_Px = dplyr::left_join(expand.grid(Intensity = vecIntensity, Season = vec.season.char, aggLevel = aggLevels), 
                         params_px, by = "Season") %>%
   mutate(Px = get_Px_Intensity_aggLevel(Intensity,aggLevel/1280,a_mu,b_mu,a_sigma,b_sigma))
  
  # Px obs  
  ggplot_px_intensity_time = ggplot(data = params_MRC$px_aggLevel_Intensity)+
    geom_point(aes(x = Intensity, y = value, color = factor(aggLevel)))+
    geom_line(data = predict_Px, aes(x = Intensity, y = Px, color = factor(aggLevel)))+
    scale_y_continuous(limits = c(0,1))+ 
    scale_x_continuous(trans = 'log10', limits = range(vecIntensity))+
    scale_color_manual(values = get_col_aggLevel(aggLevels),
                       breaks = aggLevels,
                       name = "T.scale [min]") +
    labs(x = latex2exp::TeX(r'(Precipitation Intensity [mm h$^{-1}$])'),
                   y = latex2exp::TeX(r'($P_x$)'),
                   title = "Px as function of temporal aggregation level and intensity")+
    facet_wrap(~Season)+
    theme_bw()
  
  list_plots[["Px_aggLevel_intensity"]] = ggplot_px_intensity_time
  
  int_pars_px = dplyr::left_join(fit_intermittency_mu_sigma, params_px,
                          by = "Season") %>% 
   mutate(pred_mu = a_mu*log(aggLevel/1280)+b_mu,
           pred_sigma = a_sigma*log(aggLevel/1280)+b_sigma)
  
  ggplot_px_params_aggLevel = ggplot(int_pars_px)+
    geom_line(aes(x = aggLevel, y = pred_mu))+
    geom_point(aes(x = aggLevel, y = mu, color = "mu"),fill = "black",shape = 15)+
    geom_line(aes(x = aggLevel, y = pred_sigma))+
    geom_point(aes(x = aggLevel, y = Sigma,color = "sigma"), shape = 17)+
    scale_x_continuous(trans = 'log10', breaks = aggLevels)+
    scale_color_manual(name = "",breaks = c("mu","sigma"), values=c("black","black"))+
    facet_wrap(~Season) + theme_bw()+
    labs(title = "Px params as function of time", 
         x = "Temporal aggregation level", y = "Value")+
    guides(colour = guide_legend(override.aes = list(shape=c(15,17))))
  list_plots[["Px_scaling_params"]] = ggplot_px_params_aggLevel
  
  
  # ======================= Plot: Alpha model =======================
  
  # alpha as function of temporal scale 
  
  # predict alpha star values by the model 
  predict_alpha = dplyr::left_join(expand.grid(aggLevel = aggLevels, Season = vec.season.char), 
                            params_alpha, by = "Season") %>%
   mutate(alpha = get_alpha_aggLevel(alpha0, H, res_aggLevel = aggLevel/1280))
  
  # alpha 
  ggplot_alpha_ts = ggplot(data = params_MRC$alpha_aggLevel)+
    geom_line(data = predict_alpha, aes(x = aggLevel, y = alpha))+
    geom_point(aes(x = aggLevel, y = value))+
    scale_y_continuous(trans = 'log10', limits = c(0.2,12))+ 
    scale_x_continuous(trans = 'log10', breaks = aggLevels)+
    facet_wrap(~Season)+
    labs(title = "Alpha as function of temporal scale", 
         x = "Temporal aggregation level", y = latex2exp::TeX(r'($\alpha$)'))+
    theme_bw()
  
  list_plots[["Alpha_aggLevel"]] = ggplot_alpha_ts
  
  ggplot_alpha_ts 
  
  # alpha temporal scale and intensity class 
  ggplot_alpha_ts_int = ggplot()+
    geom_point(data = params_MRC$alpha_aggLevel_Intensity,
               aes(x = Intensity, y = value, color = factor(aggLevel)))+
    scale_y_continuous(trans = 'log10', limits = c(0.2,12))+ 
    scale_x_continuous(trans = 'log10', limits = range(vecIntensity))+
    scale_color_manual(values = get_col_aggLevel(aggLevels),
                       breaks = aggLevels,
                       name = "T.scale [min]") + 
    labs(x = latex2exp::TeX(r'(Precipitation Intensity [mm h$^{-1}$])'),
         latex2exp::TeX(r'($\alpha$)'))+
    facet_wrap(~Season)+
    theme_bw()
  
  ggplot_alpha_ts_int
  list_plots[["Alpha_aggLevel_intensity"]] = ggplot_alpha_ts_int
  
  # predict alpha star values by the model 
  predict_alpha_star = dplyr::left_join(expand.grid(Intensity = vecIntensity, Season = vec.season.char), 
                                 params_alpha, by = "Season") %>%
   mutate(alpha_star = alpha_star_Intensity(Intensity,c0,c1,c2))
  
  # alpha star
  ggplot_alpha_star = ggplot()+
    geom_point(data = params_MRC$alpha_star_aggLevel_Intensity,
               aes(x = Intensity, y = value, color = factor(aggLevel)))+
    geom_line(data = predict_alpha_star,aes(x = Intensity, y = alpha_star))+
    scale_y_continuous(trans = 'log10', limits = c(0.2,12))+ 
    scale_x_continuous(trans = 'log10', limits = range(vecIntensity))+
    scale_color_manual(values = get_col_aggLevel(aggLevels),
                       breaks = aggLevels,
                       name = "T.scale [min]") + 
    labs(x = latex2exp::TeX(r'(Precipitation Intensity [mm h$^{-1}$])'),
         y = latex2exp::TeX(r'($\alpha^{*}$)'))+
    facet_wrap(~Season)+
    theme_bw()
  
  ggplot_alpha_star  
  list_plots[["Alpha_star_aggLevel_intensity"]] = ggplot_alpha_star
  
  
  # ======================= Plot: Asymm model =======================
  
  predict_phi = dplyr::left_join(expand.grid(zIndex = seq(0,1,length.out = 30), Season = vec.season.char), 
                          params_asymm, by = "Season") %>%
    mutate(phi = get_phi_z(z = zIndex, nu = nu))
  
  predict_mean_weights = dplyr::left_join(expand.grid(zIndex = seq(0,1,length.out = 30), Season = vec.season.char), 
                                   params_asymm, by = "Season") %>%
    mutate(m = get_mean_z(z = zIndex, lambda = lambda))
  
  
  
  ggplot_asymm_ratio = ggplot()+
    geom_line(data = predict_phi,aes(x = zIndex, y = phi))+
    geom_point(data = params_MRC$phi_z, aes(x = zIndex, y = value))+
    facet_wrap(~Season) + theme_bw()+
    scale_x_continuous(limits = c(0,1))+
    scale_y_continuous(limits = c(0,1))+
    labs(title = "Probability asymmetry ratio", 
         x = "z", y = latex2exp::TeX(r'($\varphi = p_{"01"}/(p_{"01"}+p_{10})$)'))
  ggplot_asymm_ratio
  
  ggplot_asymm_mean = ggplot()+
    geom_line(data = predict_mean_weights,aes(x = zIndex, y = m))+
    geom_point(data = params_MRC$emp_mean_z ,aes(x = zIndex, y = value))+
    facet_wrap(~Season) + theme_bw()+
    scale_x_continuous(limits = c(0,1))+
    scale_y_continuous(limits = c(0,1))+
    labs(title = "Mean of the cascade weights", 
         x = "z", y = "m")
  ggplot_asymm_mean
  
  list_plots[["Asymm_ratio"]] = ggplot_asymm_ratio
  list_plots[["Asymm_mean"]] = ggplot_asymm_mean
  
  params_all = as.data.frame(dplyr::left_join(params_asymm, 
                                              dplyr::left_join(params_alpha,params_px,by = "Season"),
                                       by = "Season"))
  
  rownames(params_all) = params_all$Season
  
  return(list(data = params_MRC, params = params_all, fig_plots = list_plots))
}



#==============================================================================
#' get_scaling_params_Intensity
#' 
#' Get scaling parameters of a MRC model, where the cascade weights are 
#' considered to be dependent on the temporal aggregation level 
#' and on the intensity of precipitation that has to be disaggregated.
#' In the first step observed weights are calculated for each temporal aggregation level. 
#' Also at this stage, for each observed precipitation intensity, its index z is calculated.  
#' Second, MRC parameters are estimated by considering a dependency to 
#' temporal aggregation level, precipitation intensity class and asymmetery index class z. 
#' In the last step, scaling model are fitted on estimated MRC parameter by 
#' the method of linear least square. Plots showing estimted params and their scaling models 
#' are also returned   
#' 
#' @param vecPrecip a vector of observed precipitations 
#' @param vecDates a vector of dates
#' @param resVecPrecip resolution of time series in minutes
#' @param aggLevels aggregation levels in the cascade estimation procedure, in minutes 
#' @param by_season should the parameter estimation must be done on seasonal basis 
#' @param threshold for weight discarding in all estimation procedure, in mm 
#' @param threshold_alpha_int for discarding alpha estimated at small intensities, in mm
#' @param threshold_asymm for discarding weights of small intensities for estimation of asymmetry parameters, in mm
#' 
#' 
#' @return to_do a list of parameters, plots for fitting, and df of empirical MRC params
#' 
#' @author Kaltrina Maloku
#' 
#' @export
get_scaling_params_Intensity = function(vecPrecip,
                                        vecDates,
                                        resVecPrecip = 10,
                                        by_season = T,
                                        threshold = 0,
                                        threshold_alpha_int = 0.8,
                                        threshold_asymm = 0.8,
                                        aggLevels = c(80,160,320,640,1280)){
  
  # estimate parameters of MRC 
  params_MRC = estimate_params_MRC(vecPrecip = vecPrecip, vecDates = vecDates,
                                   resVecPrecip = resVecPrecip,
                                   by_season = T,
                                   threshold = threshold, 
                                   threshold_alpha_int = threshold_alpha_int,
                                   threshold_asymm = threshold_asymm, 
                                   aggLevels = aggLevels)
  
  # some initial parameters for fitting models 
  mu_start = -0.5; sigma_start = 1
  a_start = -1;b_start = -1
  nu_start = 0.5; lambda_start = 0.5
  K_start = 0.1
  
  vec.season.char = c("DJF", "MAM", "JJA", "SON")
  
  # where to predict values 
  vecIntensity = exp(seq(from = log(0.003), to = log(150), length.out = 100))
  
  # ======================= Px model =======================
  
  # First stage: estimate params of mu and sigma for each temporal scales 
  fit_intermittency_mu_sigma = params_MRC$px_aggLevel_Intensity %>% group_by(Season) %>% 
    dplyr::do(fit = robustbase::nlrob(value~get_Px_Intensity(Intensity = Intensity, 
                                                 mu = mu, Sigma = Sigma),
                                      control = list(maxiter = 1000, printEval = F, minFactor = 1/100000),
                          start = list(mu = mu_start, Sigma = sigma_start), maxit = 2000, data = .data)) %>%
   mutate(b_mu = fit$coefficients[1], b_sigma = fit$coefficients[2],
           a_mu = 0,a_sigma = 0) %>% dplyr::select(-fit)
  
  
  # prepare results 
  params_px = fit_intermittency_mu_sigma
  
  
  
  # ======================= Alpha model ======================= 
  # estimate scaling params of mu over temporal scales 
  fit_alpha_intensity_K = params_MRC$alpha_aggLevel_Intensity %>% group_by(Season) %>%
    dplyr::do(fit = robustbase::nlrob(value~get_alpha_intensity_1par(Intensity = Intensity,
                                                         I_min = 0.1,I_max = 10, K = K),
                                      control = list(maxiter = 1000, printEval = F, minFactor = 1/100000),
                          start = list(K = K_start), 
                          data = .data, maxit = 2000)) %>% 
   mutate(K = fit$coefficients[1]) %>% dplyr::select(-fit)
  
  
  params_alpha = fit_alpha_intensity_K
  
  
  
  # ======================= Asymmetry model =======================
  # fit linear model on estimated means for each z index class
  fit_mean_weights = params_MRC$emp_mean_z %>% dplyr::filter(!is.na(zIndex)) %>% group_by(Season) %>% 
    dplyr::do(fit = robustbase::nlrob(value~get_mean_z(z = zIndex, lambda = lambda), 
                                      control = list(maxiter = 1000, printEval = F, minFactor = 1/100000),
                          start = list(lambda = lambda_start), maxit = 2000, data = .data)) %>%
   mutate(lambda = fit$coefficients[1]) %>% dplyr::select(-fit)
  
  # fit erf model on estimated ratios for each z index class
  fit_intermittency_ratio_phi = params_MRC$phi_z %>% dplyr::filter(!is.na(zIndex)) %>% group_by(Season) %>% 
    dplyr::do(fit = robustbase::nlrob(value~get_phi_z(z = zIndex, nu = nu),
                                      control = list(maxiter = 1000, printEval = F, minFactor = 1/100000),
                          start = list(nu = nu_start), maxit = 2000, data = .data)) %>%
   mutate(nu = fit$coefficients[1]) %>% dplyr::select(-fit)
  
  params_asymm = dplyr::left_join(fit_mean_weights,fit_intermittency_ratio_phi, by = "Season")
  
  # ======================= Plot: Px model =======================
  list_plots = list()
  # predict values based on the parameters of the model
  predict_Px = dplyr::left_join(expand.grid(Intensity = vecIntensity, Season = vec.season.char, aggLevel = aggLevels), 
                         params_px, by = "Season") %>%
   mutate(Px = get_Px_Intensity_aggLevel(Intensity,aggLevel/1280,a_mu,b_mu,a_sigma,b_sigma))
  
  # Px obs  
  ggplot_px_intensity_time = ggplot(data = params_MRC$px_aggLevel_Intensity)+
    geom_line(data = predict_Px, aes(x = Intensity, y = Px))+
    geom_point(aes(x = Intensity, y = value, color = factor(aggLevel)))+
    scale_y_continuous(limits = c(0,1))+ 
    scale_x_continuous(trans = 'log10', limits = range(vecIntensity))+
    scale_color_manual(values = get_col_aggLevel(aggLevels),
                       breaks = aggLevels,
                       name = "T.scale [min]") + 
    labs(x = latex2exp::TeX(r'(Precipitation Intensity [mm h$^{-1}$])'),
                  y = latex2exp::TeX(r'($P_x$)'),
                  title = "Px as function of intensity")+
    facet_wrap(~Season)+
    theme_bw()
  
  ggplot_px_intensity_time 
  
  list_plots[["Px_intensity"]] = ggplot_px_intensity_time
  
  # ======================= Plot: Alpha model =======================
  # predict alpha as function of intensity  
  predict_alpha_int = dplyr::left_join(expand.grid(Intensity = vecIntensity, Season = vec.season.char), 
                                fit_alpha_intensity_K, by = "Season") %>%
   mutate(alpha = get_alpha_intensity_1par(Intensity = Intensity,
                                            I_min = 0.1,I_max = 10, K = K))
  
  ggplot_alpha_int = ggplot()+
    geom_line(data = predict_alpha_int,aes(x = Intensity, y = alpha))+
    geom_point(data = params_MRC$alpha_aggLevel_Intensity, 
               aes(x = Intensity, y = value, color = factor(aggLevel)))+
    scale_y_continuous(trans = 'log10', limits = c(0.2,12))+ 
    scale_x_continuous(trans = 'log10', limits = range(vecIntensity))+
    scale_color_manual(values = get_col_aggLevel(aggLevels),
                       breaks = aggLevels,
                       name = "T.scale [min]") +
    labs(x = latex2exp::TeX(r'(Precipitation Intensity [mm h$^{-1}$])'),
                  y = latex2exp::TeX(r'($\alpha$)'),
                  title = "Alpha as function of intensity")+
    facet_wrap(~Season)+
    theme_bw()
  
  list_plots[["Alpha_intensity"]] = ggplot_alpha_int
  
  # ======================= Plot: Asymm model =======================
  
  predict_phi = dplyr::left_join(expand.grid(zIndex = seq(0,1,length.out = 30), Season = vec.season.char), 
                          params_asymm, by = "Season") %>%
   mutate(phi = get_phi_z(z = zIndex, nu = nu))
  
  predict_mean_weights = dplyr::left_join(expand.grid(zIndex = seq(0,1,length.out = 30), Season = vec.season.char), 
                                   params_asymm, by = "Season") %>%
   mutate(m = get_mean_z(z = zIndex, lambda = lambda))
  
  
  
  ggplot_asymm_ratio = ggplot()+
    geom_line(data = predict_phi,aes(x = zIndex, y = phi))+
    geom_point(data = params_MRC$phi_z, aes(x = zIndex, y = value))+
    facet_wrap(~Season) + theme_bw()+
    scale_x_continuous(limits = c(0,1))+
    scale_y_continuous(limits = c(0,1))+
    labs(title = "Probability asymmetry ratio", 
         x = "z", y = latex2exp::TeX(r'($\varphi = p_{"01"}/(p_{"01"}+p_{10})$)'))

  ggplot_asymm_mean = ggplot()+
    geom_line(data = predict_mean_weights,aes(x = zIndex, y = m))+
    geom_point(data = params_MRC$emp_mean_z ,aes(x = zIndex, y = value))+
    facet_wrap(~Season) + theme_bw()+
    scale_x_continuous(limits = c(0,1))+
    scale_y_continuous(limits = c(0,1))+
    labs(title = "Mean of the cascade weights", 
         x = "z", y = "m")

  list_plots[["Asymm_ratio"]] = ggplot_asymm_ratio
  list_plots[["Asymm_mean"]] = ggplot_asymm_mean
  
  params_all = as.data.frame(dplyr::left_join(params_asymm, 
                                              dplyr::left_join(params_alpha,params_px,by = "Season"),
                                       by = "Season"))
  
  rownames(params_all) = params_all$Season
  
  return(list(data = params_MRC, params = params_all, fig_plots = list_plots))
}

#==============================================================================
#' disaggregate_precip_MRC_Intensity_aggLevel
#' 
#' Disaggregate coarse resolution time series to fine resolution time series following 
#' the a MRC model where the dependency to aggregation levels and intensity is accounted for  
#' 
#' @param vecPrecip_target a vector of observed precipitations at coarser resolution
#' @param vecDates_target a vector of dates 
#' @param params_scaling a matrix of data frame of the needed parameters for disaggregation
#' @param res_coarse_aggLevel resolution of time series that needs to be disaggregated in minutes
#' @param res_fine_aggLevel target resolution of the disaggregated time series 
#' @param nb_scenarios number of disaggregated scenarios
#' @param asymmetry_option should the disaggregation be dependent on asymmetry model. Default FALSE
#' 
#' @return a matrix, each column one disaggregated scenario
#' 
#' @author Kaltrina Maloku
#' 
#' @export
disaggregate_precip_MRC_Intensity_aggLevel = function(vecPrecip_target,
                                                      vecDates_target,
                                                      params_scaling,
                                                      res_coarse_aggLevel = 1280,
                                                      res_fine_aggLevel = 40,
                                                      nb_scenarios = 10,
                                                      asymmetry_option = F){
  
  
  if(length(vecPrecip_target) != length(vecDates_target))
    stop("Vector of dates should have the same length of vector of precipitations")
  vecSeason_aggLevel = month2season(lubridate::month(vecDates_target))
  vecSeason_aggLevel = factor(vecSeason_aggLevel, levels = 1:4, 
                              labels = c("DJF","MAM","JJA","SON"))
  
  vecSeason_aggLevel = as.character(vecSeason_aggLevel)

  # a data frame of disaggregated scenarios
  data_precip = NULL
  
  #=================== Model A ==================
  set.seed(3)
  
  for (i_scenario in 1:nb_scenarios){
    message("Scenario ", i_scenario)
    i_aggLevel = res_coarse_aggLevel
    while(i_aggLevel > res_fine_aggLevel){
      # get vector of precip and seasons at the current aggregation level
      if(i_aggLevel == res_coarse_aggLevel){ # if the first level take target days 
        vecPrecip_aggLevel_current = vecPrecip_target
        vecSeason_aggLevel_current = vecSeason_aggLevel
       }
      # length of current time series 
      n_current = length(vecPrecip_aggLevel_current)
      
      # prepare data frame 
      # convert precip amounts to precipitaion intensity mm h  
      data_current = data.frame(Intensity = vecPrecip_aggLevel_current*(60/i_aggLevel),
                                zIndex = get_z_index(vecPrecip_aggLevel_current),
                                Season = vecSeason_aggLevel_current)
      
      # add params to intensity and season, each season its respective parameters
      data_current_params = dplyr::left_join(data_current, params_scaling,
                                             by = "Season")
      
      # from intesities and scaling parameters, get MRC parameters 
      # rnd_unif is random vector used for disaggregation 
      data_current_params = data_current_params %>%
       mutate(Px = get_Px_Intensity_aggLevel(Intensity,i_aggLevel/1280, a_mu,b_mu,a_sigma,b_sigma),
               alpha = get_alpha_aggLevel(alpha0, H, res_aggLevel = i_aggLevel/1280)*
                 alpha_star_Intensity(Intensity,c0,c1,c2),
               phi = get_phi_z(z = zIndex, nu = nu),
               m = get_mean_z(z = zIndex, lambda = lambda),
               rnd_unif = runif(nrow(data_current_params)),bdc = NA)
      
      
      if(asymmetry_option){    # if asymmetry option is true 
        # first determine 0/1 weights 
        data_current_params = data_current_params %>% 
         mutate(P01 = phi*(1-Px), # proba that all precip goes in the second half
                 P10 = (1-phi)*(1-Px)) %>% # proba that all precip goes in the first half
         mutate(zeros = 1*!(P01 - rnd_unif >= 0), # confort the proba to random devaites 
                 ones = 1*(1-P10 - rnd_unif < 0)) %>% # and decide if 0, 1, or in between
         mutate(bdc = replace(bdc, zeros == 0,0)) %>%
         mutate(bdc = replace(bdc, ones == 1,1)) %>% 
          dplyr::select(-c("zeros","ones","rnd_unif"))
        
        # from alpha, and the mean as modelled by asymmetry model
        # find alpha1 and alpha2
        data_current_params = data_current_params %>% 
         mutate(alpha1 = get_alpha1(alpha, m),
                 alpha2 = get_alpha2(alpha, m))
        
        indices_weights = is.na(data_current_params$bdc) & 
          !is.na(data_current_params$zIndex)
        
        # draw random deviates of beta distribution only where the weight is not yet defined
        postive_weights = mapply(function(x,y) stats::rbeta(1, shape1 = x, shape2 = y),
                                 data_current_params$alpha1[indices_weights],
                                 data_current_params$alpha2[indices_weights])
      } else {  # if asymmetry option is not true 
        # first determine 0/1 weights 
        data_current_params = data_current_params %>% 
         mutate(zeros = 1*!(0.5*(1-Px) - rnd_unif >= 0),
                 ones = 1*(0.5*(1+Px) - rnd_unif < 0)) %>% 
         mutate(bdc = replace(bdc, zeros == 0,0)) %>%
         mutate(bdc = replace(bdc, ones == 1,1)) %>% 
          dplyr::select(-c("zeros","ones","rnd_unif"))
        indices_weights = is.na(data_current_params$bdc) & 
          !is.na(data_current_params$zIndex)
        # draw random deviates of beta distribution only where the weight is not yest defined
        postive_weights = sapply(data_current_params$alpha[indices_weights],
                                 function(x) rbeta(1, shape1 = x, shape2 = x))
        
      }
      
      # replace missing weights 
      data_current_params = data_current_params %>% 
       mutate(bdc = replace(bdc, is.na(bdc) & !is.na(zIndex), postive_weights)) %>% 
       mutate(bdc = replace(bdc, is.na(bdc),0))
      
      
      # prepare a new vector of precipiation with resolution twice higher than 
      vecPrecip_aggLevel_previous = rep(NA, 2*n_current)
      
      # the first half
      vecPrecip_aggLevel_previous[seq(1, by = 2, length.out = n_current)] = 
        vecPrecip_aggLevel_current*data_current_params$bdc
      
      # to the second half 
      vecPrecip_aggLevel_previous[seq(2, by = 2, length.out = n_current)] = 
        vecPrecip_aggLevel_current*(1-data_current_params$bdc)
      
      # update aggregation level
      i_aggLevel = i_aggLevel/2
      
      # update precipitation and season vectors
      vecPrecip_aggLevel_current = vecPrecip_aggLevel_previous
      vecSeason_aggLevel_current = rep(vecSeason_aggLevel_current, each = 2)
      
    }
    
    # add to the matrix of scenarios
    data_precip = cbind(data_precip, vecPrecip_aggLevel_current)
    colnames(data_precip)[i_scenario] = paste0("result.",i_scenario)
  }
  
  # return simulated scenarios
  return(data_precip)
}

#==============================================================================
#' disaggregate_precip_MRC_Intensity
#' 
#' Disaggregate coarse resolution time series to fine resolution time series following 
#' the a MRC model where the dependency to intensity is accounted for  
#' 
#' @param vecPrecip_target a vector of observed precipitations at coarser resolution
#' @param vecDates_target a vector of dates 
#' @param params_scaling a matrix of data frame of the needed parameters for disaggregation
#' @param res_coarse_aggLevel resolution of time series that needs to be disaggregated in minutes
#' @param res_fine_aggLevel target resolution of the disaggregated time series 
#' @param nb_scenarios number of disaggregated scenarios
#' @param asymmetry_option should the disaggregation be dependent on asymmetry model. Default FALSE
#' 
#' @return a matrix, each column one disaggregated scenario
#' 
#' @author Kaltrina Maloku
#' 
#' @export
disaggregate_precip_MRC_Intensity = function(vecPrecip_target,
                                             vecDates_target,
                                             params_scaling,
                                             res_coarse_aggLevel = 1280,
                                             res_fine_aggLevel = 40,
                                             nb_scenarios = 10,
                                             asymmetry_option = F){
  
  
  if(length(vecPrecip_target) != length(vecDates_target))
    stop("Vector of dates should have the same length of vector of precipitations")
  vecSeason_aggLevel = month2season(lubridate::month(vecDates_target))
  vecSeason_aggLevel = factor(vecSeason_aggLevel, levels = 1:4, 
                              labels = c("DJF","MAM","JJA","SON"))
  
  
  # a data frame of disaggregated scenarios
  data_precip = NULL
  
  #=================== Model B ==================
  set.seed(3)
  
  for (i_scenario in 1:nb_scenarios){
    message("Scenario ", i_scenario)
    i_aggLevel = res_coarse_aggLevel
    while(i_aggLevel > res_fine_aggLevel){
      
      #message(i_aggLevel)
      
      # get vector of precip and seasons at the current aggregation level
      if(i_aggLevel == res_coarse_aggLevel){ # if the first level take target days 
        vecPrecip_aggLevel_current = vecPrecip_target
        vecSeason_aggLevel_current = vecSeason_aggLevel
      }
      # length of current time series 
      n_current = length(vecPrecip_aggLevel_current)
      
      # prepare data frame 
      # convert precip amounts to precipitaion intensity mm h  
      data_current = data.frame(Intensity = vecPrecip_aggLevel_current*(60/i_aggLevel),
                                zIndex = get_z_index(vecPrecip_aggLevel_current),
                                Season = vecSeason_aggLevel_current)
      
      # add params to intensity and season, each season its respective parameters
      data_current_params = dplyr::left_join(data_current, params_scaling,
                                             by = "Season")
      
      # from intesities and scaling parameters, get MRC parameters 
      # rnd_unif is random vector used for disaggregation 
      data_current_params = data_current_params %>%
       mutate(Px = get_Px_Intensity_aggLevel(Intensity,i_aggLevel/1280, a_mu,b_mu,a_sigma,b_sigma),
               alpha = get_alpha_intensity_1par(Intensity = Intensity,
                                                I_min = 0.1,I_max = 10, K = K),
               phi = get_phi_z(z = zIndex, nu = nu),
               m = get_mean_z(z = zIndex, lambda = lambda),
               rnd_unif = runif(nrow(data_current_params)),bdc = NA)
      
      
      if(asymmetry_option){    # if asymmetry option is true 
        # first determine 0/1 weights 
        data_current_params = data_current_params %>% 
         mutate(P01 = phi*(1-Px), # proba that all precip goes in the second half
                 P10 = (1-phi)*(1-Px)) %>% # proba that all precip goes in the first half
         mutate(zeros = 1*!(P01 - rnd_unif >= 0), # confort the proba to random devaites 
                 ones = 1*(1-P10 - rnd_unif < 0)) %>% # and decide if 0, 1, or in between
         mutate(bdc = replace(bdc, zeros == 0,0)) %>%
         mutate(bdc = replace(bdc, ones == 1,1)) %>% 
          dplyr::select(-c("zeros","ones","rnd_unif"))
        
        # from alpha, and the mean as modelled by asymmetry model
        # find alpha1 and alpha2
        data_current_params = data_current_params %>% 
         mutate(alpha1 = get_alpha1(alpha, m),
                 alpha2 = get_alpha2(alpha, m))
        
        indices_weights = is.na(data_current_params$bdc) & 
          !is.na(data_current_params$zIndex)
        
        # draw random deviates of beta distribution only where the weight is not yet defined
        postive_weights = mapply(function(x,y) stats::rbeta(1, shape1 = x, shape2 = y),
                                 data_current_params$alpha1[indices_weights],
                                 data_current_params$alpha2[indices_weights])
      } else {  # if asymmetry option is not true 
        # first determine 0/1 weights 
        data_current_params = data_current_params %>% 
         mutate(zeros = 1*!(0.5*(1-Px) - rnd_unif >= 0),
                 ones = 1*(0.5*(1+Px) - rnd_unif < 0)) %>% 
         mutate(bdc = replace(bdc, zeros == 0,0)) %>%
         mutate(bdc = replace(bdc, ones == 1,1)) %>% 
          dplyr::select(-c("zeros","ones","rnd_unif"))
        indices_weights = is.na(data_current_params$bdc) & 
          !is.na(data_current_params$zIndex)
        # draw random deviates of beta distribution only where the weight is not yet defined
        postive_weights = sapply(data_current_params$alpha[indices_weights],
                                 function(x) stats::rbeta(1, shape1 = x, shape2 = x))
        
      }
      
      # replace missing weights 
      data_current_params = data_current_params %>% 
       mutate(bdc = replace(bdc, is.na(bdc) & !is.na(zIndex), postive_weights)) %>% 
       mutate(bdc = replace(bdc, is.na(bdc),0))
      
      
      # prepare a new vector of precipitation with resolution twice higher than 
      vecPrecip_aggLevel_previous = rep(NA, 2*n_current)
      
      # the first half
      vecPrecip_aggLevel_previous[seq(1, by = 2, length.out = n_current)] = 
        vecPrecip_aggLevel_current*data_current_params$bdc
      
      # to the second half 
      vecPrecip_aggLevel_previous[seq(2, by = 2, length.out = n_current)] = 
        vecPrecip_aggLevel_current*(1-data_current_params$bdc)
      
      # update aggregation level
      i_aggLevel = i_aggLevel/2
      
      # update precipitation and season vectors
      vecPrecip_aggLevel_current = vecPrecip_aggLevel_previous
      vecSeason_aggLevel_current = rep(vecSeason_aggLevel_current, each = 2)
      
    }
    
    # add to the matrix of scenarios
    data_precip = cbind(data_precip, vecPrecip_aggLevel_current)
    colnames(data_precip)[i_scenario] = paste0("result.",i_scenario)
  }
  
  # return simulated scenarios
  return(data_precip)
}

