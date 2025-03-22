# Trying out simEventData2

func1 <- function(N) rbinom(N, 1, 0.2)
func2 <- function(N) rnorm(N)

add_cov <- list(func1, func2)
beta <- matrix(rnorm(6*4), ncol = 4, nrow = 6)

data <- simEventData2(N = 10000, add_cov = add_cov, beta = beta)
t_data <- IntFormatData(data)

# Checking that the additional covariates affect the intensities in the right way
survfit_death <- coxph(Surv(tstart, tstop, Delta == 1) ~ I(L0/50) + A0 + L + A + L1 + L2,
                       data = t_data)

survfit_oper <- coxph(Surv(tstart, tstop, Delta == 2) ~ I(L0/50) + A0 + L + L1 + L2,,
                      data = t_data[A == 0])

# benchmarking

my_bench <- bench::press(
  N = 2^(4:10),
  {
    bench::mark(
      "simEventData2" = simEventData2(N = N),
      #"simEventData2 add_cov" = simEventData2(N = N, add_cov = add_cov),
      "simEventData" = simEventData(N = N),
      check = FALSE
    )
  }
)

my_bench |>
  dplyr::mutate(expr = as.character(expression), median = as.numeric(median)) |>
  ggplot(aes(N, median, color = expr)) + geom_point() +
  #scale_y_log10() +
  geom_line() + labs(y = "time (ms)")

# Not too slow


