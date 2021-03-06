% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plot_hsmm_fit_status}
\alias{plot_hsmm_fit_status}
\title{Visualizes the status of the EM-procedure.}
\usage{
plot_hsmm_fit_status(fit_output, title = NULL, y_axis_limits = NULL)
}
\arguments{
\item{fit_output}{the output of the \code{fit_hsmm()} function.}

\item{title}{(optional) an character specifying the title of the visualization (typically the name of the model). By default, there is no title.}

\item{y_axis_limits}{(optional) a 2-element vector specifying the limits of the y-axis.}
}
\value{
a ggplot object.
}
\description{
If a model has been fitted using the function \code{fit_hsmm()}, this function can be used to visualize the convergence of the EM.
}
\examples{
my_model = simple_model
X_sim = simulate_hsmm(model = my_model, n_state_transitions = 20)
fit_results = fit_hsmm(model = my_model, X = X_sim)
plot_hsmm_fit_status(fit_output = fit_results)
}
