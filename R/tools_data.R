#' World Health Organization TB data
#'
#'A fictive 10 minute resolution time series of precipitation 
#'
#' @format ## `data_test_MRC`
#' A data frame of two columns and 1577952 rows:
#' \describe{
#'   \item{obs}{Precipitation amounts}
#'   \item{date}{Dates}
#' }
"data_test_MRC"


#' @importFrom stats runif rbeta
#' @importFrom dplyr group_by mutate summarise_at select rename filter
#' @importFrom dplyr left_join
#' @importFrom ggplot2 aes labs guides guide_legend ggplot geom_line geom_point scale_x_continuous scale_y_continuous scale_color_manual facet_wrap theme_bw
#' @importFrom magrittr %>% 
#' 
#' @export
NULL



