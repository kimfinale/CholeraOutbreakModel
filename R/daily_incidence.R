daily_incidence <- function(pars){
  # use the global variable PARAMETERS, which holds all the parameters
  # and initial conditions
  # pars only stores estimated parameters
  params <- PARAMETERS
  for (n in names(pars)) {
    params[[n]] <- pars[[n]]
  }

  out <- params$model(params)
  day_filter <- seq(1, by=round(1/params$tau), length.out=(params$ndays+1))
  mv <- params$measure_var
  out <- out[day_filter, mv, drop=FALSE]

  # if measured variables is more than one
  df <- data.frame(matrix(NA, nrow=(nrow(out)-1), ncol=length(mv)))
  names(df) <- mv
  for(i in 1:length(mv)){
    df[,i] <- diff(out[, mv[i]])
  }

  return(df)
}
