% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hsmm_functions_simulate.R
\name{simulate_hsmm}
\alias{simulate_hsmm}
\title{Simulate a state and observation sequence from a hidden semi-Markov model.}
\usage{
simulate_hsmm(
  model,
  seq_id = NULL,
  n_state_transitions = NULL,
  n_timepoints = NULL,
  seed = NULL,
  all_states = FALSE,
  min_tp_per_state = NULL,
  mult_seq_allowed = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{model}{a object of class or \code{hsmm}.}

\item{seq_id}{(optional) a \code{character} specifying the name of the sequence to simulate.}

\item{n_state_transitions}{(optional) an \code{integer} specifying the number of state transitions that should be simulated.}

\item{n_timepoints}{(optional) an \code{integer} specifying the number of time-point to simulate.
One can either specify \code{n_state_transitions} or \code{n_timepoints}.
If both are specified, \code{n_timepoints} is ignored and a sequence with \code{n_state_transitions} is specified.}

\item{seed}{(optional) the seed for the random generator. Allows to reproduce a simulated sequence.}

\item{all_states}{(optional) a \code{logical} specifying if all states should be present in the simulated sequence.
By default, this option is \code{FALSE}. When \code{TRUE}, the sequence may be longer than the specified length via \code{n_state_transitions} or \code{n_timepoints}}

\item{min_tp_per_state}{(optional) if \code{all_states == TRUE}, this option can be used to specify the minimum number of time-point in each state.}

\item{mult_seq_allowed}{(optional) if \code{all_states == TRUE}, this option can be used to specify if representation of all states in the simulated sequence can be achieved by simulated several sequences or by simulating a single sequence.}

\item{verbose}{(optional) whether to print additional info when running.}
}
\value{
A \code{data.frame} with the following columns: \code{seq_id}, the identifier of the sequence, \code{t} the time-stamp, \code{state} the hidden state and one additional column for each specified variable.
}
\description{
This function returns a \code{data.frame} with a simulated state sequence and observations.
The length of the simulated sequence can be specified either as a number of state transition (with argument \code{n_state_transitions}) or as a number of time-point (with argument \code{n_timepoints})
Each row of the returned \code{data.frame} is a time-point and the columns specify the sequence id, the time-point, the state and the values of the observation. Each variable has its own column.
}
\examples{
my_model = simple_model # simple_model is a model attached to the HiddenSemiMarkov package for demos
Xsim = simulate_hsmm(model = my_model, n_state_transitions = 10)
plot_hsmm_seq(model = my_model, X = Xsim)
}
\references{
J. O’Connell, S. Hojsgaard, Hidden semi-Markov models for multiple observation sequences: The mhsmm package for R. Journal of Statistical Software. 39, 1–22 (2011) \url{https://www.jstatsoft.org/article/view/v039i04}
}
\keyword{HSMM}
