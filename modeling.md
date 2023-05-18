Cholera outbreak model
================
2023-02-27

## Load the package

``` r
devtools::load_all()
```

## Run the model

``` r
devtools::load_all()
init <- initialize_params()

print_params <- function(params){
  n <- names(params)
  for(i in seq_along(n)){
    cat(paste0(n[i], "=", params[n[i]]), ", ")
  }
}
print_params(init)
```

``` r
PARAMETERS <- initialize_params()
# PARAMETERS$model <- seiarw
pop <- 1e5
pop_scale <- 0.07
prop_init_R <- 0.6
prop_init_I <- 0.001 
PARAMETERS$init$R <- prop_init_R * pop * pop_scale
PARAMETERS$init$I <- prop_init_I * pop * pop_scale
PARAMETERS$init$S <- pop * pop_scale - PARAMETERS$init$I - PARAMETERS$init$I

plot(1:nrow(d), d$CI, type="l")

# ?data.table::frollsum
df <- data.frame(matrix(NA, nrow=52, ncol=10))
df$week=1:52
for(i in 1:10){
  parm <- c(R0=2, R0W=1, case_threshold=10, alpha=(i-1)/10)
  d <- daily_incidence(pars=parm)
  weekroll <- data.table::frollsum(d$CI, n=7)
  sq <- seq(7, length(weekroll), by=7)
  df[, i] <- weekroll[sq]
}

head(df, n=20)
names(df) <- c(paste0(seq(0,90,10)," %"), "Week")
df <- df[,c(seq(1,9,2),11)]
library(tidyr)
# dflong <- pivot_longer(df, cols=1:ncol(df), names_to = "Reduction")
dflong <- pivot_longer(df, cols=(!Week), names_to = "Reduction")
library(ggplot2)
mycolors <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
                "black", "gold1", "skyblue2", "palegreen2", "#FDBF6F",
                "gray70", "maroon", "orchid1", "darkturquoise",
                "darkorange4", "brown")

p <- ggplot(dflong) +
  geom_line(aes(Week, value, color=Reduction), linewidth=1)+
  labs(x="Week", y="Case")+
  scale_color_manual("Intervention", values=mycolors[1:10])+
  theme_classic()+
  theme(panel.grid.major = element_line(colour="grey95", linewidth=0.3),
    text = element_text(size=14), axis.text = element_text(size=14))
p

ggsave("plots/outbreak_intervention.png", p, width=3.4*2, height=2.7*2, unit="in")



# sum(d$CI)
cat("outbreak size:", round(sum(weekly)), "\n")
# df <- data.frame(week = 1:length(weekly), case = weekly)
# library(ggplot2)
# p <- ggplot(df) +
#   geom_col(aes(week, case), fill = "brown") + 
#   labs(x="Week", y="Case")+
#   theme_classic()+
#   theme(panel.grid.major = element_line(colour="grey95", linewidth=0.3),
#     text = element_text(size=14),
#         axis.text = element_text(size=14))
# p
# ggsave("plots/typical_outbreak.png", p, width=3.4*2, height=2.7*2, unit="in")

out <- PARAMETERS$model(PARAMETERS)
day_filter <- seq(1, by=round(1/PARAMETERS$tau), length.out=(PARAMETERS$ndays+1))
out <- out[day_filter, ]
```

## Vaccination campaign

``` r
PARAMETERS <- initialize_params()
# PARAMETERS$model <- seiarw
pop <- 1e5
pop_scale <- 0.07
prop_init_R <- 0.6
prop_init_I <- 0.001 
PARAMETERS$init$R <- prop_init_R * pop * pop_scale
PARAMETERS$init$I <- prop_init_I * pop * pop_scale
PARAMETERS$init$S <- pop * pop_scale - PARAMETERS$init$I - PARAMETERS$init$I

parm <- c(R0=2, R0W=1, case_threshold=10, alpha=0, vacc_cov=0.4, 
          vacc_eff_v1 = 0.9, vacc_eff_v2 = 0.9)
params <- PARAMETERS
for (n in names(parm)) {
  params[[n]] <- parm[[n]]
}

out <- seiarw(params)

d <- daily_incidence(pars=parm, variable = names(out))
plot(1:nrow(d), d$CI, type="l")
```

### Plot

``` r
# ?data.table::frollsum
df <- data.frame(matrix(NA, nrow=52, ncol=10))
df$week=1:52
for(i in 1:10){
  parm <- c(R0=2, R0W=1, case_threshold=10, alpha=(i-1)/10)
  d <- daily_incidence(pars=parm)
  weekroll <- data.table::frollsum(d$CI, n=7)
  sq <- seq(7, length(weekroll), by=7)
  df[, i] <- weekroll[sq]
}

head(df, n=20)
names(df) <- c(paste0(seq(0,90,10)," %"), "Week")
df <- df[,c(seq(1,9,2),11)]
library(tidyr)
# dflong <- pivot_longer(df, cols=1:ncol(df), names_to = "Reduction")
dflong <- pivot_longer(df, cols=(!Week), names_to = "Reduction")
library(ggplot2)
mycolors <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
                "black", "gold1", "skyblue2", "palegreen2", "#FDBF6F",
                "gray70", "maroon", "orchid1", "darkturquoise",
                "darkorange4", "brown")

p <- ggplot(dflong) +
  geom_line(aes(Week, value, color=Reduction), linewidth=1)+
  labs(x="Week", y="Case")+
  scale_color_manual("Intervention", values=mycolors[1:10])+
  theme_classic()+
  theme(panel.grid.major = element_line(colour="grey95", linewidth=0.3),
    text = element_text(size=14), axis.text = element_text(size=14))
p

ggsave("plots/outbreak_intervention.png", p, width=3.4*2, height=2.7*2, unit="in")



# sum(d$CI)
cat("outbreak size:", round(sum(weekly)), "\n")
# df <- data.frame(week = 1:length(weekly), case = weekly)
# library(ggplot2)
# p <- ggplot(df) +
#   geom_col(aes(week, case), fill = "brown") + 
#   labs(x="Week", y="Case")+
#   theme_classic()+
#   theme(panel.grid.major = element_line(colour="grey95", linewidth=0.3),
#     text = element_text(size=14),
#         axis.text = element_text(size=14))
# p
# ggsave("plots/typical_outbreak.png", p, width=3.4*2, height=2.7*2, unit="in")

out <- PARAMETERS$model(PARAMETERS)
day_filter <- seq(1, by=round(1/PARAMETERS$tau), length.out=(PARAMETERS$ndays+1))
out <- out[day_filter, ]
```

### SEIARW Erlang distributed model

## 

surveillance activity?

``` r
# extract outbreak characteristics (outbreak duration, peak incidence)
out1=out
out2=out+rnorm(nrow(out))

nsamp = 1
nfeat = 3
threshold = 7
popsize = 1e5
df = data.frame(matrix(NA, nrow=nsamp, ncol=nfeat))
names(df) = c("peak_inc","outbreak_dur","attack_rate")
for(i in 1:nsamp){
  weekly = diff(out$CI.1)  
  df[i, "peak_inc"] = max(out$CI.1)
  over = which(weekly > threshold)
  df[i, "outbreak_dur"] = out$time[max(over)] - out$time[min(over)]
  df[i, "attack_rate"] = tail(out$CI.1, 1) / popsize
}
df
```

# Grid search

See if how the proportion of already-immune people for a given threshold

``` r
params <- initialize_params()
params$vacc_cov <- c(0.0, 0.0)
params$delay_until_2nd_campaign <- 30
params$alpha <- 0
params$tau <- 0.1
params$R0 = 2
params$R0W = 2

outbk = simulate_outbreak(params=params)

run_outbreak = function(params, 
                        model=seiarw_2ag_erlang,
                        popsize=1e5,
                        inc_ref_pop=1e3,
                        prop_s=1, 
                        prop_eff_pop=1,
                        output_var=c("CI.1","CI.2"),
                        output_time=7,
                        threshold=1){
  
  params$init$R <- (1-prop_s) * popsize * prop_eff_pop
  params$init$I <- 0.001 * popsize * prop_eff_pop
  params$init$S <- popsize * prop_eff_pop - params$init$I

  out <- run_model(model, params, 
                   output_time=output_time,
                   output_var=output_var)
  
  df = data.frame(ci = out$CI.1 + out$CI.2)
  df$time = out$time / output_time # make it a week
  
  weekly = diff(df$ci)
  over = which(weekly > threshold) # index for the week of which the incid is over the threshold
  start = df$time[min(over)]
  end = df$time[max(over)] 
  peak = df$time[which.max(weekly)]
  peak_inc_per_1e3 = weekly[peak] / popsize * inc_ref_pop
  time_to_peak = (peak - start)
  outbreak_dur = (end - start)
  outbreak_size = tail(df$ci, 1)   
  attack_rate_per_1e3 = outbreak_size / popsize * inc_ref_pop   
  
  return(c(outbreak_dur=outbreak_dur, 
           time_to_peak=time_to_peak,
           peak_inc_per_1e3=peak_inc_per_1e3, 
           outbreak_size=outbreak_size,
           attack_rate_per_1e3=attack_rate_per_1e3))
}

res = run_outbreak(params=params)

get_outbreak_char = function(df, threshold=1, popsize=1e5,
                                        inc_ref_pop=1e3){
  weekly = diff(df$ci)
  over = which(weekly > threshold)
  start = df$time[min(over)]
  end = df$time[max(over)] 
  peak = df$time[which.max(weekly)]
  peak_inc_per_1e3 = weekly[peak] / popsize * inc_ref_pop
  time_to_peak = (peak - start)
  outbreak_dur = (end - start)
  outbreak_size = tail(df$ci, 1)   
  attack_rate_per_1e3 = outbreak_size / popsize * inc_ref_pop   
  
  return(c(outbreak_dur=outbreak_dur, time_to_peak=time_to_peak,
              peak_inc_per_1e3=peak_inc_per_1e3, 
              outbreak_size = outbreak_size,
              attack_rate_per_1e3=attack_rate_per_1e3))
}

out <- run_model_with_varying_S(prop_s=0.3)
plot(out$time, c(0,diff(out$ci)), type="b")
get_outbreak_char(df=out)


plot(out$time, c(0,diff(out$ci)), type="b")
cnt <- 1
my_discrete_colors[1:9]
for(i in seq(0.1,0.9,0.1)){
  
  out <- run_model_with_varying_S(params=params,prop_s=1-i)
  lines(out$time, c(0,diff(out$ci)), col=my_discrete_colors[cnt], lwd=2)
 
  res = get_outbreak_char(df=out, threshold=7, popsize=1e5)
  cat("dur=", res$outbreak_dur, ", size=", res$outbreak_size, ", time_to_peak=", res$time_to_peak, ", peak_inc=", 
      res$peak_inc_per_1e3, ", ar=", res$attack_rate_per_1e3, "\n")
  cnt = cnt + 1
}
```

### EpiEstim

``` r
params <- initialize_params()
params$vacc_cov <- c(0.0, 0.0)
params$delay_until_2nd_campaign <- 30
params$alpha <- 0
params$tau <- 0.01
params$R0W = 0
params$prop_children <- 0
params$fA = 0
params$sigma = 0 # or immunity waning
params$R0=2

out=simulate_outbreak(params=params,
                      prop_s=1,
                      output_time=1)

library(EpiEstim)
day_incid = diff(out$timeseries$ci)
t_start <- seq(2, length(day_incid)-1)   
t_end <- t_start + 1                 
res <- estimate_R(incid = day_incid, 
                  method = "parametric_si",
                  config = make_config(list(
                  mean_si = 4.5,
                  std_si = 3, 
                  t_start = t_start, 
                  t_end = t_end)))
res$R$`Mean(R)`[1:20]

plot(res$R$`Mean(R)`[1:50])
```

## Grid search

``` r
params <- initialize_params()
params$vacc_cov <- c(0.0, 0.0)
params$delay_until_2nd_campaign <- 30
params$alpha <- 0
params$tau <- 0.1
params$R0 <- 2
params$R0W <- 0
# first few columns are parameters varied
res = expand.grid(prop_s=seq(0.1,1.0,0.1), 
                  threshold=1,
                  eff_pop_prop=seq(0.1,0.9,0.1))

res$outbreak_dur = NA
res$outbreak_size = NA
res$time_to_peak = NA
res$peak_inc_per_1e3 = NA
res$attack_rate_per_1e3 = NA
res$prop_peak_inc = NA

output = vector("list", nrow(res))

for(i in 1:nrow(res)){
  cat("i =", i, "\n")
  out <- simulate_outbreak(params=params, 
                           prop_s=res$prop_s[i],
                           prop_eff_pop=res$eff_pop_prop[i])
  
  output[[i]] = out
  res[i, "outbreak_dur"] = out$outbreak_char[["outbreak_dur"]]
  res[i, "outbreak_size"] = out$outbreak_char[["outbreak_size"]]
  res[i, "time_to_peak"] = out$outbreak_char[["time_to_peak"]]
  res[i, "peak_inc_per_1e3"] = out$outbreak_char[["peak_inc_per_1e3"]]
  res[i, "attack_rate_per_1e3"] = out$outbreak_char[["attack_rate_per_1e3"]] 
  res[i, "prop_peak_inc"] = out$outbreak_char[["prop_peak_inc"]] 

}
```

### $R_t = R_0 S$ during the initial period

$R_t$ measured over the initial period of an outbreak through the use of
EpiEstim package gives $R_0 S$ where $S$ is the fraction of susceptibles
rather than $R_0$. For an endemic disease like cholera, $S$ can be lower
than 1 and therefore, $R_t$ during the initial phase of a cholera
outbreak can lead to an underestimate of $R_0$. See the code block below
to check if the incidence is the same regardless of $R_0$ as long as the
$R_0 S$ remains constant. This implies that given the value of initial
Rt we can choose different combinations of $R_0$ and $S$. Instead, we
assumed that $R_t=R_0$ and assumed $S=1$. Then, we explored the varying
size of the population (i.e., S) that produced the observed
characteristics of the outbreak.

``` r
R0S = 4
parm <- initialize_params()
parm$vacc_cov <- c(0.0, 0.0)
parm$delay_until_2nd_campaign <- 30
parm$alpha <- 0
parm$tau <- 0.1
parm$R0 <- 2
parm$R0W <- 0
prop_s <- R0S/parm$R0 

out = simulate_outbreak(params=parm, prop_s=prop_s, output_time=1)
head(out$timeseries, n=10)


parm2 <- initialize_params()
parm2$vacc_cov <- c(0.0, 0.0)
parm2$delay_until_2nd_campaign <- 30
parm2$alpha <- 0
parm2$tau <- 0.1
parm2$R0 <- 4
parm2$R0W <- 0
prop_s = R0S/parm2$R0

out2 = simulate_outbreak(params=parm2, prop_s=prop_s, output_time=1)
head(out2$timeseries, n=20)

df=data.frame(ci1=out$timeseries$ci, ci2=out$timeseries$ci)
df$diff <- df$ci1 - df$ci2
summary(df$diff)
```
