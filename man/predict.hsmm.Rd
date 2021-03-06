% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hsmm_functions.R
\name{predict.hsmm}
\alias{predict.hsmm}
\title{Predicts the hidden state sequence from observations}
\usage{
\method{predict}{hsmm}(object, newdata, method = "Viterbi", verbose = FALSE, ...)
}
\arguments{
\item{object}{an \code{hsmm} model specified via the \code{specify_hsmm()} function.}

\item{newdata}{a \code{data.frame} with the observation sequences.}

\item{method}{a \code{character} specifying the method to be used, i.e. either \code{"Viterbi"} or \code{"Fwbw"}.}

\item{verbose}{a logical (default = \code{FALSE}). Should the function print additional information?}

\item{...}{additional arguments.}
}
\description{
This function is a wrapper around the function \code{predict_states_hsmm()}.
It predicts the most likely hidden states from observations.
Two methods are implemented:
\code{"Viterbi"} applies the Viterbi algorithm and predicts the most likely sequence of hidden states,
and \code{"FwBw"} applies the Forward-Backward algorithm and returns the probability of each state at each time-point.
}
\seealso{
see \code{predict_states_hsmm()} for the full description, options and examples.
}
