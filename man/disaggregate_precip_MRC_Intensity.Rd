% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools_scaling_MRC.R
\name{disaggregate_precip_MRC_Intensity}
\alias{disaggregate_precip_MRC_Intensity}
\title{disaggregate_precip_MRC_Intensity}
\usage{
disaggregate_precip_MRC_Intensity(
  vecPrecip_target,
  vecDates_target,
  params_scaling,
  res_coarse_aggLevel = 1280,
  res_fine_aggLevel = 40,
  nb_scenarios = 10,
  asymmetry_option = F
)
}
\arguments{
\item{vecPrecip_target}{a vector of observed precipitations at coarser resolution}

\item{vecDates_target}{a vector of dates}

\item{params_scaling}{a matrix of data frame of the needed parameters for disaggregation}

\item{res_coarse_aggLevel}{resolution of time series that needs to be disaggregated in minutes}

\item{res_fine_aggLevel}{target resolution of the disaggregated time series}

\item{nb_scenarios}{number of disaggregated scenarios}

\item{asymmetry_option}{should the disaggregation be dependent on asymmetry model. Default FALSE}
}
\value{
a matrix, each column one disaggregated scenario
}
\description{
Disaggregate coarse resolution time series to fine resolution time series following 
the a MRC model where the dependency to intensity is accounted for
}
\author{
Kaltrina Maloku
}
