% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forecasting_practical.R
\name{fit_seir}
\alias{fit_seir}
\title{fit the SEIR model to data}
\usage{
fit_seir(incidence_data, current_week, starting_week, R_0_min, R_0_max)
}
\arguments{
\item{incidence_data}{data frame extracted by \code{extract_incidence}}

\item{current_week}{numeric vector of length 1: week number of the current week}

\item{starting_week}{numeric vector of length 1:
guess for week number when the epidemic started. Use data from starting week
to current week to predict}

\item{R_0_min}{numeric vector of length 1: lower bound of R_0 values over
which to search}

\item{R_0_max}{numeric vector of length 1: upper bound of R_0 values over
which to search}
}
\value{
the value of R_0 with the maximum likelihood
}
\description{
\code{fit_seir} fits an SEIR model to incidence data
}
