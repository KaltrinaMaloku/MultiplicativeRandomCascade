###===============================###===============================###
### Kaltrina Maloku
### 23/02/2023, Grenoble
### IGE
### kaltrina.maloku@univ-grenoble-alpes.fr
###
### Provide tools to estimate the parameters of MRC models 
###
### add refs, rupp and my paper 
###===============================###===============================###




#==============================================================================
# return a vector of seasons 
#==============================================================================
get.list.season = function(){
  return(c('DJF','MAM','JJA','SON'))
}

#==============================================================================
#' month2season
#'
#' transform vector of months to seasons
#'
#' @param vecMonth a vector of months given as integers 1:12
#'
#' @author Guillaume Evin
month2season = function(vecMonth){
  iSeason = c(1,1,2,2,2,3,3,3,4,4,4,1)
  return(iSeason[vecMonth])
}

#==============================================================================
#' get_intensity_class_id
#' 
#' Get intensity id class given the precip amount
#' 
#' @param x a vector of precipitation amounts 
#' @param I_min Minimal observed precipitation in mm 
#' @param res_aggLevel current aggregation level in minutes 
#' 
#' @return a vector of the same length as x indicating the class of each precip amount in mm/h  
#' 
#' @author Kaltrina Maloku
#' 
#' @references Rupp, D. E., Keim, R. F., Ossiander, M., Brugnach, M., and Selker, J. S.:
#' Time scale and intensity dependency in multiplicative cascades for temporal rainfall disaggregation
#' Water Resources Research, 45,https://doi.org/10.1029/2008WR007321, 2009.
#' 
#' @export
get_intensity_class_id = function(x, I_min = 0.1,res_aggLevel){
  # if precip amount is 0 make NA 
  x[x < I_min] = NA
  
  # find n such that I_min*2^(n-1) <= x < I_min*2^n
  n_int_class = floor(log(x/I_min)/log(2))+1
  
  # the maximal classes of intensities 
  max_int_class = 15
  
  # at hourly resolution we have 
  int_bounds_intervals = I_min*2^(1:max_int_class-1)
  int_class_id = rep(NA, max_int_class-1)
  for (i in 1:(length(int_bounds_intervals)-1)){
    int_class_id[i] = sqrt(int_bounds_intervals[i]*int_bounds_intervals[i+1])
  }
  
  # at the aggregation level of interest we have 
  int_class_id_current_aggLevel = int_class_id*(60/res_aggLevel)
  
  # find the id with the help of factor 
  int_id = as.numeric(as.character(factor(n_int_class, 
                                          levels = 1:(max_int_class-1),
                                          labels = int_class_id_current_aggLevel)))
  
  return(int_id)
}

#==============================================================================
#' get_z_index
#' 
#' Get class of a given precip amount
#' 
#' @param x a vector of precipitation amounts
#' 
#' @return a vector of the same length as x, indicating asymmetry index z. 
#' It is calculated as \deqn{z_t =\frac{x_{t-1}+0.5 x_t}{x_{t-1}+x_t+x_{t+1}}}
#' 
#' @author Kaltrina Maloku
#' 
#' @references My_paper
#' 
#' @export
get_z_index = function(x){
  n.x = length(x)
  
  # ratios and the bins 
  R_2 = x[2:(n.x-1)]
  R_1 = x[1:(n.x-2)]
  R_3 = x[3:n.x]
  
  Z = (0.5+R_1/R_2)/(1+R_1/R_2+R_3/R_2)
  Z = c(0.5, Z, 0.5)
  
  # if the nearby values are unknown assign 0.5
  # this arrives in the cse os missing data 
  Z[!is.na(x) & is.na(Z)] = 0.5
  Z[x == 0] = NA
  
  return(Z)
}

#==============================================================================
#' get_class_z_bin
#' 
#' For a given z index, return the index of its bin class
#' 
#' @param z a vector of z indices  
#' @param bins_nb number of bins 
#' 
#' @return a numeric vector of the same length as x, indicating z index class 
#' 
#' @author Kaltrina Maloku
#' 
#' @references My_paper
#' 
#' @export
get_class_z_bin = function(z,bins_nb = 10){
  
  # create ten bins [0,0.1),[0.1,0.2),...,[0.9,1)
  bins_intervals = seq(0,1, length.out = bins_nb+1)
  
  # the center of the bin for each bin class
  bins_index = seq(from = 0.5*bins_intervals[2], by = bins_intervals[2] - bins_intervals[1],
                   length.out = bins_nb)
  
  # find at which bin each z index falls
  x_bins = findInterval(z,bins_intervals,left.open = T,all.inside = T)
  
  # find the bin and return it as a character 
  x_bins_factor = factor(x_bins, levels = 1:bins_nb, labels = bins_index)
  x_bins_id = as.numeric(as.character(x_bins_factor))
  
  return(x_bins_id)
}

#==============================================================================
#' sum_fixed_window
#' 
#' Sum over a fixed window  
#' 
#' @param x a vector precipitations 
#' @param k the length of fixed window `sum`
#' 
#' @return vector of length \code{floor(length(x)/k)}, sum over the window of length \code{k} 
#' @details Sums are assigned to the most left element of the window
#' @author Kaltrina Maloku
#' 
#' @export
sum_fixed_window = function(x,k){
  # check if the vector and window size are correct
  if(!class(x) %in% c("numeric","integer")) stop("x must be a numeric vector")
  if(length(x) < k) stop("length of fixed window must be smaller than the vector")
  
  # fixed sums  
  x_sum = RcppRoll::roll_sum(x = x,n = k, by = k,
                             fill = NULL, 
                             align = "left", na.rm = F) 
  
  return(x_sum)
}



#==============================================================================
#' fit_symm_beta_dist
#' 
#' Fit symmetric beta distribution, with probability density function
#' \deqn{f(x;\alpha) = \frac{1}{B(\alpha, \alpha)} {x}^{(\alpha-1)}(1-x)^{(\alpha-1)}} 
#' 
#' @param x a vector of weights 
#' 
#' @return estimated \eqn{\alpha} parameter
#' @details If the sample size is smaller than 10, return NA 
#' @author Kaltrina Maloku
#' 
#' @export
fit_symm_beta_dist = function(x){
  sample_length = length(x)
  # if less than 10 elements on the set return NA
  # do not estimate alpha
  if(sample_length < 10){
    alpha_par = NA
  } else {
    # calculate variance 
    var_pw = stats::var(x, na.rm = T)
    
    # if it happens that all data have the same value return NA
    var_pw[var_pw == 0] = NA
    
    # estimate alpha 
    alpha_par = 1/(8*var_pw)-0.5
    
    
    # if sample size is smaller than 20 perform a test to check the quality of fit
    if(sample_length < 20 & sample_length > 10 & !is.na(alpha_par)){
      probabilities = (1:sample_length)/(sample_length+1)
      pred = stats::qbeta(probabilities,
                          shape1 = alpha_par, shape2 = alpha_par)
      k = stats::ks.test(x,pred)
      if(k$p.value <= 0.1){
        print(k); print(x)
        alpha_par = NA
      }
    }
  }
  
  return(alpha_par)
}


#==============================================================================
#' get_alpha1
#' 
#' A function that estimates the parameter alpha1 of beta distribution 
#' given the parameter alpha of symmetric beta distribution and its mean 
#' 
#' @param alpha a numeric or a  numeric vector of the parameter alpha of symmetric beta distribution
#' @param m a numeric or a  numeric vector of the mean of the distribution
#' @return a numeric of the same length as \code{length(alpha)} with the parameter alpha1 of beta distribution
#' @author Kaltrina Maloku
#' 
#' @export
get_alpha1 = function(alpha, m){
  var_weights = 1/(8*alpha+4)
  mean_weights = m
  shape1 = mean_weights*((mean_weights*(1-mean_weights)/var_weights)-1)
  #shape2 = (1-mean_weights)*((mean_weights*(1-mean_weights)/var_weights)-1)
  return(shape1)
}

#==============================================================================
#' get_alpha2
#' 
#' A function that estimates the parameter alpha2 of beta distribution 
#' given the parameter alpha of symmetric beta distribution and its mean 
#' 
#' @param alpha a numeric or a  numeric vector of the parameter alpha of symmetric beta distribution
#' @param m a numeric or a  numeric vector of the mean of the distribution
#' @return a numeric of the same length as \code{length(alpha)} with the parameter alpha2 of beta distribution
#' @author Kaltrina Maloku
#' 
#' @export
get_alpha2 = function(alpha, m){
  var_weights = 1/(8*alpha+4)
  mean_weights = m
  #shape1 = mean_weights*((mean_weights*(1-mean_weights)/var_weights)-1)
  shape2 = (1-mean_weights)*((mean_weights*(1-mean_weights)/var_weights)-1)
  return(shape2)
}

#==============================================================================
#' get_portion_of_0
#' 
#' Finds portion of cascade weights that are equal to 0 
#' 
#' @param x a vector of weights 
#' @return portion of weights equal to 0 
#' @details If the sample size is smaller than 10, return NA 
#' @author Kaltrina Maloku
#' 
#' @export
get_portion_of_0 = function(x){
  if(length(x)<10){
    pr_0 = NA
  } else {
    x = x[!is.na(x)]
    pr_0 = mean(x == 0)
  }
  return(pr_0)
}

#==============================================================================
#' get_portion_of_1
#' 
#' Finds portion of cascade weights that are equal to 1
#' 
#' @param x a vector of weights 
#' @return portion of weights equal to 1
#' @details If the sample size is smaller than 10, return NA 
#' @author Kaltrina Maloku
#' 
#' @export
get_portion_of_1 = function(x){
  if(length(x)<10){
    pr_1 = NA
  } else {
    x = x[!is.na(x)]
    pr_1 = mean(x == 1)
  }
  return(pr_1)
}

#==============================================================================
#' get_portion_of_1
#' 
#' Finds portion of cascade weights that are strictly greater than 0 and strictly smaller than 1
#' 
#' @param x a vector of weights 
#' @return portion of weights that are strictly greater than 0 and strictly smaller than 1 
#' @details If the sample size is smaller than 10, return NA 
#' @author Kaltrina Maloku
#' @export
get_portion_of_x = function(x){
  if(length(x)<10){
    pr_x = NA
  } else {
    x = x[!is.na(x)]
    pr_x = mean(x < 1 & x > 0)
  }
  return(pr_x)
}

#==============================================================================
#' get_Px_Intensity
#' 
#' Intermittency model, it returns the probability to divide precipitation intensity
#' given the parameters of the model \eqn{\mu} and \eqn{\sigma} 
#' 
#' @param Intensity a numeric value, the precipitation intensity of interest 
#' @param mu,Sigma parameters of the model
#' 
#' @return a numeric value 
#' 
#' @author Kaltrina Maloku
#' 
#' @references Rupp, D. E., Keim, R. F., Ossiander, M., Brugnach, M., and Selker, J. S.:
#' Time scale and intensity dependency in multiplicative cascades for temporal rainfall disaggregation
#' Water Resources Research, 45,https://doi.org/10.1029/2008WR007321, 2009.
#'
#' @export
get_Px_Intensity = function(Intensity, mu, Sigma){
  E = pracma::erf((log(Intensity) - mu)/(sqrt(2)*Sigma))
  return(0.5+0.5*E)
}


#==============================================================================
#' get_Px_Intensity_aggLevel
#' 
#' Intermittency model, it returns the probability to divide precipitation intensity
#' given the parameters of the model \eqn{a_\mu,b_\mu,a_\sigma,b_\sigma} and 
#' current aggregation level
#' 
#' @param Intensity a numeric vector or value, the precipitation intensity of interest 
#' @param a_mu,b_mu,a_sigma,b_sigma parameters of the model
#' @param res_aggLevel a numeric value, current aggregation level in minutes
#' 
#' @return a numeric value or a vector the same length as \code{length(Intensity)}
#' 
#' @author Kaltrina Maloku
#' 
#' @references Rupp, D. E., Keim, R. F., Ossiander, M., Brugnach, M., and Selker, J. S.:
#' Time scale and intensity dependency in multiplicative cascades for temporal rainfall disaggregation
#' Water Resources Research, 45,https://doi.org/10.1029/2008WR007321, 2009.
#'
#' @export
get_Px_Intensity_aggLevel = function(Intensity,res_aggLevel,a_mu,b_mu,a_sigma,b_sigma){
  mu = a_mu*log(res_aggLevel)+b_mu
  Sigma = a_sigma*log(res_aggLevel)+b_sigma
  E = pracma::erf((log(Intensity) - mu)/(sqrt(2)*Sigma))
  return(0.5+0.5*E)
}

#==============================================================================
#' alpha_star_Intensity
#' 
#' The quadratic model of the scaled \eqn{\alpha}. Given the parameters of the model \code{c0, c1, c2} and the 
#' precipitation intensity, it returns the value of scaled \code{alpha}. 
#' 
#' @param Intensity a numeric vector or value, the precipitation intensity of interest 
#' @param c0,c1,c2 the parameters of the model
#' 
#' @return a numeric value or a vector the same length as \code{length(Intensity)}
#' 
#' @author Kaltrina Maloku
#' 
#' @references Rupp, D. E., Keim, R. F., Ossiander, M., Brugnach, M., and Selker, J. S.:
#' Time scale and intensity dependency in multiplicative cascades for temporal rainfall disaggregation
#' Water Resources Research, 45,https://doi.org/10.1029/2008WR007321, 2009.
#'
#' @export
alpha_star_Intensity = function(Intensity, c0, c1, c2){
  return(exp(c0+c1*log(Intensity)+c2*log(Intensity)*log(Intensity)))
}


#==============================================================================
#' get_alpha_aggLevel
#' 
#' Model of alpha as function of temporal scale It returns the parameter \eqn{\alpha} of beta distribution 
#' given the parameters of the model \eqn{\alpha_{0}} and \eqn{H}
#' 
#' @param res_aggLevel a numeric value, current aggregation level in minutes
#' @param alpha0,H parameters of the model
#' 
#' @return a numeric value, the value of the parameter \eqn{\alpha} given aggregation level  
#' 
#' @author Kaltrina Maloku
#' 
#' @references Rupp, D. E., Keim, R. F., Ossiander, M., Brugnach, M., and Selker, J. S.:
#' Time scale and intensity dependency in multiplicative cascades for temporal rainfall disaggregation
#' Water Resources Research, 45,https://doi.org/10.1029/2008WR007321, 2009.
#'
#' @export
get_alpha_aggLevel = function(alpha0, H, res_aggLevel){
  return(alpha0*res_aggLevel^(H))
}

#==============================================================================
#' get_alpha_intensity_1par
#' 
#' A model of the parameter \eqn{\alpha} as function of the precipitation intensity.
#' Its only parameter explains the degree of dependency to the intensity. 
#' The model is considered to be constant for intensities smaller than \code{I_min}, 
#' and higher than \code{I_max}.  
#' 
#' @param Intensity a numeric vector or value, the precipitation intensity of interest 
#' @param I_min for intensities smaller than this value the model is constant. Default value: 0.1 mm/h 
#' @param I_max for intensities higher than this value the model is constant. Default value: 10 mm/h 
#' @param K the parameter of the model 
#' 
#' @return a value or numeric vector of length \code{length(Intensity)} value 
#' 
#' @author Kaltrina Maloku
#' 
#' @references my_paper
#' @export
get_alpha_intensity_1par = function(Intensity, I_min, I_max, K){
  I = log(Intensity)
  I1 = log(I_min)
  I2 = log(I_max)
  
  alpha = (1*(I<I1)) + (exp(K*(I-I1)^2))*(I>=I1 & I<=I2) + (exp(K*(I2-I1)^2))*(I>I2)
  
  return(alpha)
}


#==============================================================================
#' get_phi_z
#' 
#' Probability asymmetry ratio model, it returns the ratio
#' \eqn{\varphi = \frac{p_{01}{p_{01}+p{10}}}} given the precipitation asymmetry indez \code{z} 
#' and the parameter \code{nu}
#' 
#' @param z a vector or a numeric value of z index  
#' @param nu parameter of the model
#' 
#' @return a value or numeric vector of length \code{length(z)}
#' 
#' @author Kaltrina Maloku
#' 
#' @references my_paper
#' @export
get_phi_z = function(z, nu){
  E = pracma::erf((0.5-z)/(sqrt(2)*nu))
  return(0.5+0.5*E)
}

#==============================================================================
#' get_mean_z
#' 
#' Model of the mean of the cascade weights as function of asymmetry index z. 
#' Given the parameter of the model, \code{lambda}, and the asymmetry index z it returns 
#' modelled mean of the distribution. 
#' 
#' @param z a vector or a numeric value of z index  
#' @param lambda parameter of the model
#' 
#' @return a value or numeric vector of length \code{length(z)} value 
#' 
#' @author Kaltrina Maloku
#' 
#' @references my_paper
#' @export
get_mean_z = function(z, lambda){
  return(lambda*(z-0.5)+0.5)
}



#==============================================================================
#' get_col_aggLevel
#' 
#' Intermittency model, it returns the parameter \eqn{\alpha} of beta distribution 
#' given the parameters of the model \eqn{\alpha_{0}} and \eqn{H}
#' 
#' @param aggLevel a vector providing the aggregation levels in minutes
#' 
#' @return a vector of colors, same length as \code{length(aggLevel)}
#' 
#' @author Kaltrina Maloku
#' 
#' @export
get_col_aggLevel = function(aggLevel){
  col_aggLevel = viridis::viridis(length(aggLevel), direction = -1)
  return(col_aggLevel)
}

#==============================================================================
#' get_indices_center_1280
#' 
#' A function that helps extract 1280 minute "days" 
#' from a time series of 1440 minute days.  
#' Given the number of days, it returns the indices in 10 min data of a centered day at 12:00
#' under the assumption that the days are consecutive
#' 
#' @param nb_days a numeric, number of days in the time series 
#' @return a numeric vector of indices  
#' @author Kaltrina Maloku
#' 
#' @examples
#' get_indices_center_1280(nb_days = 2)
#' 
#' @export
get_indices_center_1280 = function(nb_days){
  i = 1:nb_days
  ind_1280 = as.vector(sapply(i, function(x){
    x = (144*(x-1)+8+1):(144*(x)-8)
  }))
  return(ind_1280)
}

#==============================================================================
#' estimate_params_MRC
#' 
#' Get parameters of MRC as estimated on observed weights
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
#' @return to_do 
#' 
#' 
#' 
#' @export
estimate_params_MRC = function(vecPrecip,vecDates,
                               resVecPrecip = 10,
                               by_season = T,
                               threshold,
                               threshold_alpha_int = 0.8,
                               threshold_asymm = 0.8,
                               aggLevels = c(40,80,160,320,640,1280)){
  
  #  ======================== Preliminary checks  ========================
  # preliminary verifications of vector types and lengths 
  if(length(vecPrecip) != length(vecDates)) stop("Length of dates and precipitaions must be the same")
  if(!any(c("POSIXlt","POSIXct") %in% class(vecDates))) stop("vecDates must be a vector of class date")
  if(!class(vecPrecip) %in% c("numeric","integer")) stop("vecPrecip must be a numeric vector")
  
  # Temporal resolution of time series  
  if(is.null(resVecPrecip)){# it resolution is not given find it from dates vecto
    resVecPrecip = difftime(vecDates[2], vecDates[1], units = "mins")
    resVecPrecip = as.numeric(resVecPrecip)
  } else { # verify if it corresponds to what is declared 
    if(resVecPrecip != difftime(vecDates[2],vecDates[1], units = "mins"))
      stop("Temporal resolution of time series does not correspond with what is declared")
  }
  
  # check if aggregation levels are correct 
  if(!any(round(aggLevels/(2*resVecPrecip)) == aggLevels/(2*resVecPrecip))) stop("aggLevels must be multiple of 2*resVecPrecip")
  
  # check if the time differences are correct 
  diff_time = difftime(vecDates[2:length(vecDates)], vecDates[1:(length(vecDates)-1)])
  diff_time = as.numeric(diff_time)
  
  # if continuous time series of 10-minute resolution remove first 80 min 
  # and remove last 80 mins of each day 
  if(all(diff_time == resVecPrecip)){ # if time difference for the whole time series is 10 minutes
    message("Removing 80 mins from the begining and the end of each day")
    # number of days
    nb_days = round(length(vecDates)/(24*60/resVecPrecip))
    # get indices that correspond to the   
    indices_1280min = get_indices_center_1280(nb_days)
    vecDates = vecDates[indices_1280min]
    vecPrecip = vecPrecip[indices_1280min]
  } else { # make sure that each "day" has 128 observations
    if(!all(unique(diff_time) %in% c(10, 170,NA))){ # a first filter, the difference between observations
      nb_obs_per_day = table(as.Date(vecDates))
      
      if(!all(nb_obs_per_day[-length(nb_obs_per_day)] == 128)){
        stop("At least a day does not have the right number of observations
             the number of observation on each day must be 128 for a time series
             of 10 minute resolution")
        #print(nb_obs_per_day[which(nb_obs_per_day != 128)])
      }
    }
  }
  
  if(by_season){
    vecSeason = month2season(lubridate::month(vecDates))
  } else {
    vecSeason = 1
  }
  
  #  ======================== Start aggregating  ========================
  i_aggLevel = aggLevels[1]
  DF_weights_obs = NULL
  for(i_aggLevel in aggLevels){
    # fixed sum of the previous aggregation level
    vecPrecip_aggLevel_previous = sum_fixed_window(x = vecPrecip, k = 0.5*i_aggLevel/resVecPrecip)
    
    # fixed sum of the current aggregation level
    vecPrecip_aggLevel_current = sum_fixed_window(x = vecPrecip_aggLevel_previous, k = 2)
    
    n_current = length(vecPrecip_aggLevel_current)
    
    # precip amount on the first half
    vecPrecip_aggLevel_previous_1st_half = vecPrecip_aggLevel_previous[seq(from = 1, by = 2,
                                                                           length.out = n_current)]
    
    # observed weights 
    weights_current = vecPrecip_aggLevel_previous_1st_half/vecPrecip_aggLevel_current
    
    # vector of season for the current level
    vecSeason_current = vecSeason[seq(from = 1, by = i_aggLevel/resVecPrecip, 
                                      length.out = n_current)]
    
    # find z index at the current aggrgation level 
    vec_z_index = get_z_index(vecPrecip_aggLevel_current)
    
    # find z class id 
    vec_z_index_id = get_class_z_bin(vec_z_index)
    
    # find intensity class id at the current aggregation level
    vec_intensity_id =  get_intensity_class_id(vecPrecip_aggLevel_current, I_min = 0.1, 
                                               res_aggLevel = i_aggLevel)
    
    
    # create a data frame that holds all the necessary information 
    DF_aggLevel_current = data.frame(Season = vecSeason_current,
                                     Precip = vecPrecip_aggLevel_current,
                                     Weight = weights_current,
                                     Intensity = vec_intensity_id,
                                     zIndex = vec_z_index_id,
                                     aggLevel = i_aggLevel) %>%
       filter(!is.na(Weight) & Precip > threshold) %>% dplyr::select(-"Precip")
    
    DF_weights_obs = rbind(DF_weights_obs,DF_aggLevel_current)
    
    
  }
  
  DF_weights_obs = DF_weights_obs %>% dplyr::filter(!is.na(Season))
  DF_weights_obs$Season = factor(DF_weights_obs$Season, levels = 1:4,
                                 labels = c("DJF", "MAM", "JJA", "SON"))
  # only positive weights
  DF_weights_obs_pos = DF_weights_obs %>% dplyr::filter(Weight != 0 & Weight != 1)
  
  #  ======================== Alpha parameter ========================
  # fit symmetric Beta Distribution on observed weights for each temporal aggregation Level 
  # depending on season 
  # depending on the aggregation level
  alpha_estimates_by_aggLevel = DF_weights_obs_pos %>% 
    group_by(Season, aggLevel) %>% 
    dplyr::summarise_at("Weight", function(x) fit_symm_beta_dist(x)) %>% stats::na.omit() %>% 
    mutate(param = "alpha", zIndex = NA, Intensity = NA) %>% dplyr::rename(value = Weight) %>%
    dplyr::select(Season, param, value, aggLevel, Intensity, zIndex)
  
  # fit symmetric Beta Distribution on observed weights 
  # depending on season 
  # depending on the aggregation level
  # depending on the intensity class
  alpha_estimates_by_aggLevel_Intensity = DF_weights_obs_pos %>% 
    group_by(Season, aggLevel, Intensity) %>% 
    dplyr::summarise_at("Weight", function(x) fit_symm_beta_dist(x)) %>% stats::na.omit() %>% 
    mutate(param = "alpha", zIndex = NA) %>%dplyr::rename(value = Weight) %>%
    mutate(small_int = Intensity > threshold_alpha_int*60/aggLevel) %>%
    dplyr::filter(small_int) %>% dplyr::select(-small_int) %>%
    dplyr::select(Season, param, value, aggLevel, Intensity, zIndex)
  
  # divide each alpha estimated for temporal aggregation level and intensity class 
  # by alpha estimated for temporal scale when ignoring the intensity classes
  # remove three smallest int class for each temp agg level 
  alpha_star_estimates_by_aggLevel_Intensity = dplyr::left_join(alpha_estimates_by_aggLevel_Intensity,
                                                         alpha_estimates_by_aggLevel %>% 
                                                          dplyr::rename(alpha_t = value) %>% dplyr::select(Season,aggLevel,alpha_t),
                                                         by = c("Season", "aggLevel")) %>%
    mutate(small_int = Intensity > threshold_alpha_int*60/aggLevel) %>%
    dplyr::filter(small_int) %>% dplyr::select(-small_int) %>% mutate(value = value/alpha_t) %>%
    dplyr::select(Season, param, value, aggLevel, Intensity, zIndex)
  
  #  ======================== Intermittency parameters ======================== 
  # estimate Px, proba to divide precip into both sub-steps
  # by season, intensity and temporal aggregation level
  px_estimates_by_aggLevel_Intensity = DF_weights_obs %>% 
    group_by(Season, aggLevel, Intensity) %>% 
    dplyr::summarise_at("Weight",function(x) get_portion_of_x(x)) %>% stats::na.omit() %>% 
    mutate(param = "px", zIndex = NA) %>% dplyr::rename(value = Weight) %>%
    dplyr::select(Season, param, value, aggLevel, Intensity, zIndex)
  
  #  ======================== Asymmetry parameters ======================== 
  # estimate the mean of postive weight distribution,
  # for asymmetry model of weight distrbution
  # by season and z index class
  mean_estimates_by_z_index = DF_weights_obs_pos %>% 
    mutate(small_int = Intensity > threshold_asymm*60/aggLevel) %>%
    dplyr::filter(small_int) %>% dplyr::select(-small_int) %>%
    group_by(Season, zIndex) %>% 
     summarise_at("Weight", function(x) mean(x, na.rm = T)) %>% 
     rename(value = Weight) %>%
    mutate(param = "emp_mean", aggLevel = NA, Intensity = NA) %>%
    dplyr::select(Season, param, value, aggLevel, Intensity, zIndex) 
  
  # proba that all precip is distributed to the second half
  # by season and z index class
  p01_estimates_by_z_index = DF_weights_obs %>% 
    mutate(small_int = Intensity > threshold_asymm*60/aggLevel) %>%
     filter(small_int) %>% dplyr::select(-small_int) %>% 
    group_by(Season, zIndex) %>% 
     summarise_at("Weight",function(x) get_portion_of_0(x)) %>%
     rename(p01 = Weight)
  
  # proba that all precip is distributed to the first half
  # by season and z index class
  p10_estimates_by_z_index = DF_weights_obs %>% 
    mutate(small_int = Intensity > threshold_asymm*60/aggLevel) %>%
     filter(small_int) %>% dplyr::select(-small_int) %>%
    group_by(Season, zIndex) %>% 
     summarise_at("Weight",function(x) get_portion_of_1(x)) %>%
     rename(p10 = Weight)
  
  # estimate phi, asymmetry parameter for intermittency model
  # by season and z index class
  phi_estimates_by_z_index =  left_join(p01_estimates_by_z_index, p10_estimates_by_z_index, by = c("Season", "zIndex")) %>%
    mutate(value = p01/(p01+p10), param = "phi", aggLevel = NA, Intensity = NA) %>%
     dplyr::select(Season, param, value, aggLevel, Intensity, zIndex)
  
  # list of all estimated parameters 
  list_params_MRC = list(alpha_aggLevel = alpha_estimates_by_aggLevel,
                         alpha_aggLevel_Intensity = alpha_estimates_by_aggLevel_Intensity,
                         alpha_star_aggLevel_Intensity = alpha_star_estimates_by_aggLevel_Intensity,
                         px_aggLevel_Intensity = px_estimates_by_aggLevel_Intensity,
                         emp_mean_z = mean_estimates_by_z_index,
                         phi_z = phi_estimates_by_z_index)  
  
  
  
  return(list_params_MRC)
}

