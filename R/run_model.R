run_model <- function(model, params, output_time=1, output_var=NULL){
  # use the global variable PARAMETERS, which holds all the parameters
  # and initial conditions
  # pars only stores estimated parameters
  # params <- PARAMETERS
  # for (n in names(pars)) {# update the parameters
  #   params[[n]] <- pars[[n]]
  # }
  out <- model(params)
  time_filter <- seq(1, by=ceiling(output_time/params$tau),
                     length.out=(ceiling(params$ndays/output_time)))

  if(is.null(output_var)) {
    output_var <- params$measure_var
  }
  out <- out[time_filter, output_var, drop=FALSE]

  return(out)
}
