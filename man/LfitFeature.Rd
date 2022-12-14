% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FitFeature.R
\name{LfitFeature}
\alias{LfitFeature}
\title{LASSO fit feature
use LASSO regression optimize feature}
\usage{
LfitFeature(x, y, od = NULL, time = 1, ...)
}
\arguments{
\item{x}{feature matrix}

\item{y}{observations, survival status}

\item{od}{output dirname}

\item{time}{run times
default 1}

\item{...}{other parameters}
}
\value{
choose_feature
}
\description{
LASSO fit feature
use LASSO regression optimize feature
}
