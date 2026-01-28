# SIMULATOR -----------------------------------------------------------------

## All the functions needed to produce a simulation and all the information
## about it.

#' Extract model parameters from a simulation scenario
#'
#' @param obj Scenario created with simulate_point_release()
#'
#' @return A list with model parameters l0, D = diffusivity, a = alpha and
#'  b = beta, and the reparameterisations pds, pc0 and pc5.
#' @export
#'
#' @examples
#' simulation_parameters(sit_scenario)
#'
simulation_parameters <- function(obj) {
  ans <- with(
    parms(obj),
    list(
      l0 = -log(pds),
      D = diffusivity,
      pc0 = pc0,
      pc5 = pc5,
      a = pc0,
      b = (log(pc0) - log(pc5)) / 5**2
    )
  )
  
  ## alternative parameters
  ans <- c(
    ans,
    pds = parms(obj)$pds,
    diffusivity = parms(obj)$diffusivity,
    alpha = ans$a,
    beta = ans$b
  )
  return(ans)
  
}


#' Trap coordinates
#'
#' Retrieve trap coordinates for a simulation object
#'
#' @param x Scenario created with simulate_point_release()
#'
#' @return A 3-column tibble with trap id and x-y coordinates.
#'
#' @export
#'
#' @examples
#' trap_coords(sit_scenario)
trap_coords <- function(x) {
  ans <- parms(x)$traps.loc
  
  return(ans)
}

#' Coordinates and direction of active individuals by day
#'
#' Get results from the simulation of the diffusion and capture of insects
#' after a point release.
#'
#' @param x Scenario created with simulate_point_release()
#'
#' @return Give the evolution of the position of each individual. There are also
#' two additional attributes :
#'          -  loc: the location of the traps.
#'          -  trap_counts: counts (and cumulative counts) of the capture time
#'             for each trap.
#'
#' @export
#'
#' @examples
#' sim_results(sit_scenario)
sim_results <- function(x) {
  
  sim_outcome <- out(x)
  
  traps <- parms(x)$traps.loc
  
  ntraps <- nrow(traps)
  
  index_to_time <- function(t) (as.numeric(t) - 1) * parms(x)$tstep
  
  ## Positions of active individuals by time index
  positions <- bind_rows(
    sim_outcome,
    .id = 't'
  ) |>
    mutate(t = index_to_time(t)) |>
    tibble()
  
  ## Trap counts
  trap_counts <- lapply(
    sim_outcome,
    attr,
    'trap_counts'
  ) |>
    map(~data.frame(id = traps$id, n = .)) |>
    list_rbind(names_to = "t") |>
    mutate(
      t = index_to_time(t)
    ) |>
    group_by(id) |>
    arrange(id, t) |>
    mutate(
      n_cum = n,
      n = diff(c(0, n_cum))
    ) |>
    ungroup()
  
  attr(trap_counts, 'loc') <- traps
  
  attr(positions, 'trap_counts') <- trap_counts
  
  return(positions)
  
}


#' Get trap counts from a simulation of the diffusion and capture of insects
#' after a point release
#'
#' Trap (id) counts (n) and cumulated counts (n_cum) by day (t) with trap
#' coordinates (x, y) and distance to release point (z)
#'
#' @param x Scenario created with simulate_point_release()
#'
#' @return A tibble with the number and cumulative number of catch per trap and
#' time step (as defined in the simulation)
#'
#' @export
#'
#' @examples
#' trap_counts(sit_scenario)
trap_counts <- function(x) {
  attr(sim_results(x), 'trap_counts')
}


#' Get capture times and trap for caught individuals
#'
#' @param x Scenario created with simulate_point_release()
#'
#' @return A tibble with the number and cumulative number of catch per trap and
#' time step (as defined in the simulation)
#'
#' @export
#'
#' @examples
#' capture_times(sit_scenario)
capture_times <- function(x) {
  
  catches <- trap_counts(x)
  
  
  capture <- tibble(
    id = factor(inverse.rle(
      structure(
        list(
          lengths = catches$n,
          values = catches$id
        ),
        class = 'rle'
      )
    )),
    t = inverse.rle(
      structure(
        list(
          lengths = catches$n,
          values = catches$t
        ),
        class = 'rle'
      )
    )
  )
  
  ## Assume capture time is the center of the observation interval
  ## (currently counts are associated to the initial time)
  time_steps <- sort(unique(diff(sort(unique(catches$t)))))
  
  ## These should be all essentially the same, up to numerical error
  stopifnot(diff(range(time_steps)) < 1e-11)
  time_step <- mean(time_steps)
  
  capture$t <- capture$t + time_step / 2
  
  return(capture)
  
}


#' Simulate the diffusion and capture of insects after a point release.
#'
#' @param traps Data.frame with 3 columns. Trap id, and coordinates x and y in
#'   m, relative to the release point at the origin.
#' @param n Integer > 0. Number of insects released.
#' @param diffusivity Numeric > 0. Diffusion coefficient in m²/day
#' @param pds Numeric, between 0 and 1. Probability of Daily Survival.
#' @param pc0 Numeric, between 0 and 1. Probability of catching an insect that
#'   is just over the trap during a day.
#' @param pc5 Numeric, between 0 and 1. Probability of catching an insect that
#'   is at 5m from the trap during a day.
#' @param ndays Integer > 0. Number of days to carry out the simulation.
#' @param tstep Numeric > 0. Temporal resolution of the simulation. In
#'   days.
#'
#' @return
#' @import simecol
#' @export
#'
#' @examples
simulate_point_release <- function(
    traps,
    n,
    diffusivity,
    pds,
    pc0,
    pc5,
    ndays,
    tstep
) {
  
  sim_pars <- list(
    ninds = n,
    diffusivity = diffusivity,
    pds = pds,
    traps.loc = traps,
    pc0 = pc0,
    pc5 = pc5,
    ndays = ndays,
    tstep = tstep
  )
  
  ans <- new(
    "sitModel",
    main = function(time, init, parms) {
      init <- move(init, parms)
      init <- survive(init, parms)
      init <- capture(init, parms)
      init
    },
    parms = sim_pars,
    initfunc = function (obj) {
      
      ## Release-point coordinates
      x0 <- y0 <- 0
      
      ## Number of released individuals
      ninds <- parms(obj)$ninds
      
      ## Table of initial coordinates of individuals
      indivs <- data.frame(
        x = rep(x0, ninds),
        y = rep(y0, ninds)
      )
      
      ## Trap locations and initial counts (n)
      attr(indivs, "trap_counts") <- rep(0L, nrow(parms(obj)$traps.loc))
      
      init(obj) <- indivs
      
      obj
    },
    times = c(from = 0, to = sim_pars$ndays, by = sim_pars$tstep),
    solver = "iteration", # or default method: "iteration"
    equations = list()
  )
  
  equations(ans) <- list(
    move = function(x, parms) {
      ## D = e^2 / 4 / dt for a Brownian Motion in 2D
      n <- nrow(x)
      step_size <- sqrt(2 * parms$diffusivity * parms$tstep * rexp(n, rate = 0.5))
      a <- 2 * pi * runif(n)
      dx <- step_size * cos(a)
      dy <- step_size * sin(a)
      x$x <- x$x + dx
      x$y <- x$y + dy
      return(x)
    },
    survive = function(x, parms) {
      ## A daily survival test is conducted for each individual independently
      ans <- x[runif(nrow(x)) < parms$pds**(parms$tstep), ]
      return(ans)
    },
    capture = function(x, parms) {
      
      trap_counts <- attr(x, "trap_counts")
      stopifnot(!is.null(trap_counts))
      
      if (nrow(x) > 0) {
        ## For each trap, evaluate the probability of catching each individual around
        ## using an exponential kernel of the form a*exp(-b*x^2) where x is the
        ## distance and the parameters a and b are adjusted to the given values
        aa <- parms$pc0
        bb <- (log(aa) - log(parms$pc5)) / 5**2
        
        ## Matrix of squared distances from every individual to every trap
        D2 <- outer(x$x, parms$traps.loc$x, '-') ** 2 +
          outer(x$y, parms$traps.loc$y, '-') ** 2
        
        ## Daily capture rate matrix for each individual at each trap
        p_mat <- exp(log(aa) - bb * D2)
        
        ## Probability of no-capture in any of the traps, for each individual
        p0 <- exp(rowSums(log(1 - parms$tstep * p_mat)))
        
        ## Matrix of multinomial probabilities for sampling
        p_multim <- cbind(p0, p_mat / rowSums(p_mat) * (1 - p0))
        
        ## Fix some numerical NaNs due to 0/0
        p_multim[is.nan(p_multim)] <- 0
        
        ## Simulate captures.
        ## Indicator matrix of the outcome for each individual.
        ## Matrix transposed: (1 + n_traps) x n_individuals
        ## First row corresponds to not-captured.
        captures <- apply(
          p_multim,
          1,
          \(.) stats::rmultinom(n = 1, size = 1, prob = .)
        )
        
        not_captured <- captures[1, ] > 0
        
        ans <- x[not_captured,]
        
        trap_counts <- trap_counts + rowSums(captures[-1, , drop = FALSE])
      } else {
        ans <- x
      }
      attr(ans, 'trap_counts') <- trap_counts
      
      return(ans)
    }
  )
  
  return(simecol::sim(ans))
}


# FIGURES -----------------------------------------------------------------

plot_brownian_motion <- function(n_steps = 1e4, n_days = 4, seed = 20240703) {
  brownian_motion(n_steps, n_days, seed) |>
    ggplot() +
    facet_wrap(~day, labeller = "label_both") +
    aes(x, y, colour = day) +
    geom_path(
      data = ~rename(., facet = day),
      aes(colour = facet),
      alpha = 0.2
    ) +
    geom_path() +
    coord_fixed() +
    theme_void()
}


#' Empirical and theoretical capture times
#'
#' Histograms with recorded capture times by trap from a simulation object,
#' with the corresponding theoretical probability density functions.
#'
#' @param x A sitModel object
#' @param min_catch Positive integer. Minimum number of total catch in a trap
#' in order to be displayed in the figure.
#'
plot_capture_times <- function(x, min_catch = 10) {
  
  ## Empirical results for sit_scenario
  captures_by_trap <- capture_times(x) |>
    count(id)
  
  ## Select the trap ids with non-trivial captures
  target_traps <- captures_by_trap |>
    filter(n > 10) |>
    pull(id) |>
    as.character()
  
  max_time <- ceiling(max(capture_times(x)$t))
  
  
  ## Theoretical results for sit_scenario
  loc <- trap_coords(x) |>
    select(-id) |>
    as.matrix()
  pars <- simulation_parameters(x)
  
  # Survival end specific hazard functions
  S <- overall_survival(loc, pars)
  h <- expected_specific_hazard_functions(loc, pars)
  
  time_step <- 0.1
  time_points <- seq(0, max_time, time_step)
  
  capture_densities <-
    ## Specific hazard rates
    h(time_points) |>
    as_tibble() |>
    filter(
      i %in% target_traps
    ) |>
    inner_join(
      tibble(
        t = time_points,
        s = S(time_points)
      ),
      by = 't'
    ) |>
    group_by(i) |>
    mutate(
      z = h * s,
      ## Normalize density within traps
      f = z / sum(time_step * z)
    ) |>
    ungroup()
  
  
  
  capture_times(x) |>
    filter(id %in% target_traps) |>
    ggplot() +
    facet_wrap(~ id) +
    geom_histogram(
      aes(x = t, y = after_stat(density)),
      alpha = 0.6,
      fill = "darkorange",
      bins = 15
    ) +
    geom_line(
      data = capture_densities |> rename(id = i),
      aes(x = t, y = f),
      color='steelblue'
    ) +
    labs(
      x = "Time (days)",
      y = "Density"
    ) +
    theme_ipsum(grid = "Y")
  
}


#' Total capture probabilities
#'
#' Empirical realisation and 95% CI, vs theoretical expectation
#'
#' @param obj sitModel object.
#'
#' @examples
#' plot_capture_prob(sit_scenario)
plot_capture_prob <- function(obj) {
  
  capture_prob_data(obj) |>
    ## Look at frequencies in traps that have at least some captures
    filter(freq > 0) |>
    ggplot() +
    aes(i, prob_theoretical) +
    geom_point(
      colour = "darkorange",
      position = position_nudge(x = .15)
    ) +
    # geom_line() +
    geom_pointrange(
      aes(y = freq, ymin = inf, ymax = sup),
      colour = "grey30"
    ) +
    labs(
      x = "Trap",
      y = "Log-probability"
    ) +
    scale_y_log10() +
    theme_ipsum(grid = "Y")
}

#' Total capture probabilities
#'
#' Empirical realisation vs theoretical expectation and 95% frequency
#' interval
#'
#' @param obj sitModel object.
#'
#' @examples
#' plot_capture_prob_quantile(sit_scenario)
plot_capture_prob_quantile <- function(obj) {
  
  capture_prob_data(obj) |>
    mutate(
      i = fct_recode(i,"💀" = "0")
    ) |>
    ## Look at frequencies in traps that have at least some captures
    filter(freq > 0) |>
    ggplot() +
    aes(i, freq) +
    geom_point(
      position = position_nudge(x = .15),
      colour = "darkorange"
    ) +
    # geom_line() +
    geom_pointrange(
      aes(y = prob_theoretical,
          ymin = quant_inf / 50e3,
          ymax = quant_sup / 50e3),
      colour = "steelblue"
    ) +
    labs(
      x = "End cause",
      y = "Log-probability"
    ) +
    scale_y_log10() +
    theme_ipsum(grid = "Y")
  
}


#' Empirical and theoretical densities of position
#'
#' A 2x2 density plot of the position of insects comparing the theoretical and
#' realised distributions for 2 simulated scenarios.
#'
#' @param x A sitModel object with the realisation of one simulation.
#' @param y A sitModel object with the realisation of another simulation.
#' @param t Positive number. Time (in days).
#' @param n Positive integer. Number of grid points in each direction. Can be
#'   scalar or a length-2 integer vector.
#' @param limits Numeric vector of length 4. (xmin, xmax, ymin, ymax). Calculated
#'   by default from the ranges of variation in the data.
#'
#' @return
#' @export
plot_empirical_vs_theoretical_density <- function(obj1, obj2, t, n, limits) {
  
  var_theo <- function(obj, t) 2 * parms(obj)$diffusivity * t
  
  
  bind_rows(
    `No traps` = empirical_density(obj2, t, n, limits = limits) |>
      mutate(
        z_Theoretical = dnorm(x, sd = sqrt(var_theo(obj2, t))) *
          dnorm(y, sd = sqrt(var_theo(obj2, t)))
      ),
    `With traps` = empirical_density(obj1, t, n, limits = limits) |>
      mutate(
        z_Theoretical = dnorm(x, sd = sqrt(var_theo(obj1, t))) *
          dnorm(y, sd = sqrt(var_theo(obj1, t)))
      ),
    .id = "scenario"
  ) |>
    rename(z_Empirical = z) |>
    pivot_longer(
      cols = starts_with("z_"),
      names_prefix = "z_",
      names_to = "type",
      values_to = "z"
    ) |>
    ggplot() +
    aes(x, y) +
    facet_grid(scenario ~ type) +
    geom_raster(
      aes(fill = z),
      show.legend = FALSE
    ) +
    scale_fill_viridis_c(option = "inferno") +
    coord_fixed() +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_ipsum(grid = "") +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    )
}



# ANIMATION ------------------------------------------------------------------

#' Save an animation of a diffusion process
#'
#' @param x Data frame with variables t (time in days), x, y (coordinates)
#' @param x Data frame with variables t (time in days), x, y (coordinates) and
#' n_cum (cumulative captures)
#' @param file Output file
#'
#' Makes a animated plot using ggplot and gganimate.
#'
#' @value Returns the requested path of the saved file.
sit_animate <- function(x, traps, file) {
  anim <- ggplot(x, aes(x, y)) +
    geom_point(alpha = 0.7, show.legend = FALSE, na.rm = TRUE) +
    geom_point(
      data = traps |>
        left_join(
          attr(traps, 'loc'),
          by = 'id'
        ),
      aes(size = n_cum, group = t),
      col = 'darkorange',
      shape = 15
    ) +
    coord_fixed() +
    # theme_void() +
    labs(
      title = 'Day: {round(frame_time, 1)}',
      size = "Catch"
    ) +
    transition_time(t) +
    theme_void()
  
  
  vv <- gganimate::animate(
    anim,
    nframes = n_distinct(x$t),
    detail = 1,
    fps = 10,
    renderer = gifski_renderer()
  )
  
  gganimate::anim_save(file, animation = vv)
  return(file)
}

# MODEL FUNCTIONS ------------------------------------------------------------

## All the basic function from the known model (i.e. model where we know the
## the capture and death time).

#' Overall survival function
#'
#' @param loc Numeric 2-column matrix with trap coordinates.
#' @param pars List with parameters of the system: pds, diffusivity, alpha and
#'   beta.
#' @param n_eval Integer. Approximate lower bound for the number of evaluations
#'   points for the overal hazard function for the computation of the cumulative
#'   hazards.
#'
#' @return A function that takes a vector of times and returns a numeric vector of overall survival at the corresponding points
#' @export
#'
#' @examples
#' loc <- matrix(5 * runif(2 * 5), ncol = 2)
#' pars <- list(pds = 0.8, diffusivity = 20, alpha = 0.8, beta = 0.2)
#' s <- overall_survival(loc, pars)
#' s(1:3)
#' system.time(s(seq(0, 50, length = 4e3)))
#' tibble(t = seq(0, 3, length.out = 100)) |>
#' mutate(y = s(t)) |>
#' ggplot() +
#' aes(t, y) +
#' geom_line() +
#' scale_y_continuous(limits = c(0,1))
overall_survival <- function(loc, pars, n_eval = 100) {
  
  h <- expected_total_hazard_function(loc, pars)
  
  ans <- function(t) {
    
    t0 <- sort(unique(c(0, t)))
    partitions_intervals <- seq_interpolate(t0, length.out = n_eval)
    
    h_evals <- lapply(partitions_intervals, h)
    
    exponents <- cumsum(vapply(h_evals, \(x) simpsons_rule(x[['h']], diff(x[['t']])[1]), 1))
    
    ans <- exp(-exponents)
    
    ## If 0 was the first time point, return 1 for the survival at 0
    if (length(t) == length(t0)) ans <- c(1, ans)
    
    
    return(ans)
  }
  
  return(ans)
}




#' Expected total hazard function
#'
#' @param loc Numeric 2-column matrix with trap coordinates.
#' @param pars List with parameters of the system: pds, diffusivity, alpha and
#'   beta.
#'
#' @return A function that takes a vector of times and returns a tidy tibble
#'   with total hazard rates for time point.
#' @export
#'
#' @examples
#' loc <- matrix(5 * runif(2 * 5), ncol = 2)
#' pars <- list(pds = 0.8, diffusivity = 20, alpha = 0.8, beta = 0.2)
#' h <- expected_total_hazard_function(loc, pars)
#' h(seq(0, 3, by = 0.01)) |>
#' ggplot(aes(t, h)) +
#' geom_line()
expected_total_hazard_function <- function(loc, pars) {
  
  h_spec <- expected_specific_hazard_functions(loc, pars)
  
  ans <- function(t) {
    aggregate(h_spec(t), by = h~t, 'sum')
  }
  return(ans)
}

#' Expected specific hazard functions
#'
#' Hazard rate function for capture in each of a set of traps, including the
#' hazard rate of mortality (with number 0).
#'
#' @param loc Numeric 2-column matrix with trap coordinates.
#' @param pars List with parameters of the system: pds, diffusivity, alpha and
#'   beta.
#'
#' @return A function that takes a vector of times and returns a tidy tibble
#'   with hazard rates for each termination cause and time.
#' @export
#'
#' @examples
#' loc <- matrix(5 * runif(2 * 5), ncol = 2)
#' ## Arrange the traps in ascending order of distance from the origin
#' loc <- loc[order(rowSums(loc**2)), ]
#' pars <- list(pds = 0.3, diffusivity = 5, alpha = 8, beta = 0.2)
#' h <- expected_specific_hazard_functions(loc, pars)
#' h(seq(0, 3, by = 0.01)) |>
#' ggplot() +
#' #facet_wrap(~i, scales = 'free') +
#' aes(t, h, colour = i) +
#' geom_line()
expected_specific_hazard_functions <- function(loc, pars) {
  
  ## Hazard rate for death
  h_0 = -log(pars$pds)
  
  ## Hazard rate for capture on trap i
  gamma <- 4 * pars$diffusivity * pars$beta
  
  ans <- function(t) {
    denom <- 1 + gamma * t
    
    norm_squared <- diag(loc %*% t(loc))
    
    h_i <- exp(
      log(pars$alpha) - log(denom) - t(outer(pars$beta * norm_squared, denom, FUN = '/'))
    )
    
    ans <- data.frame(
      t = rep(t, times = 1 + nrow(loc)),
      i = factor(rep(seq.int(nrow(loc) + 1) - 1, each = length(t))),
      h = as.numeric(cbind(h_0, h_i))
    )
    
    return(ans)
    
  }
  
  return(ans)
  
}


#' Expected hazard rate for a trap at a given distance
#' kappa fixed at a given value
exp_spec_fun <- function(phi, D, r, k = 12) {
  
  function(t) {
    phi / (k + 4 * D * t) * exp(- r^2 / (k + 4 * D * t))
  }
  
}

## Parameter transformations
par_alpha <- function(phi, kappa) {
  # phi / kappa
  phi / kappa
}

par_beta <- function(phi, kappa) {
  # pi / 4 / kappa**2
  1 / kappa
}

par_phi <- function(alpha, beta) {
  # sqrt(pi) * alpha / 2 / sqrt(beta)
  alpha / beta
}

par_kappa <- function(alpha, beta) {
  # sqrt(pi) / 2 / sqrt(beta)
  1 / beta
}


# LIKELIHOOD COMPUTATION -----------------------------------------------------

## All the functions used to compute the log-likelihood of the experiment (with
## cumulative densities).

#' Log-likelihood generator
#'
#' Produces a log-likelihood function from the location of traps, counts of
#' captures by trap and period and total number of releases.
#'
#' @param loc Numeric 2-column matrix with trap coordinates.
#' @param dat Data frame with capture data. See details for structure.
#' @param N Integer. Total number of released individuals.
#'
#' @return A function that takes 4 parameter values (daily death rate,
#'   diffusivity (m²/day), daily capture rate at at trap's location and decay
#'   parameter) and returns the log-likelihood of such parameter set.
#' @export
#'
#' @examples
#' l0 = -log(0.8); D = 20; a = 0.8; b = 0.2
#' loc <- matrix(5 * runif(2 * 5), ncol = 2)
#' dat <- data.frame(
#'   i = rep(1:5, times = 3),
#'   ts = rep(0:2, each = 5),
#'   te = rep(1:3, each = 5),
#'   n = ceiling(exp(rnorm(5 * 3)))
#' )
#' N = 1e4
#' ## Define the log-likelihood function for a given scenario
#' llf <- logLik_gen(loc, dat, N)
#'
#' ## Evaluate the log-likelihood of a certain parameter set
#' ## Takes several minutes
#' llf(l0, D, a, b)
#'
#' ## llf is vectorised and works in parallel using {future}
#' plan(multisession, workers = 3)
#' llf(l0 + c(-.1, 0, .1), D = 19:21, a = 7:9/10, b = 0.2 + -1:1/1e2)
logLik_gen <- function(loc, dat, N) {
  
  ll <- function(l0, D, a, b) {
    
    stopifnot(l0 > 0, D > 0, a > 0, b > 0)
    pars <- list(
      pds = exp(-l0),
      diffusivity = D,
      alpha = a,
      beta = b
    )
    
    ll0 <- logLik_non_obs(loc, pars)  ## takes a few seconds
    li_f <- logLik_obs(loc, pars)
    
    ## Vector of log-likelihoods for each single observation in the data
    lli <- li_f(dat[[1]], dat[[2]], dat[[3]])  ## takes a few seconds
    
    ## Log-likelihood of the corresponding number of observations
    llis <- lli * dat[[4]]
    
    ## Total log-likelihood
    ans <- ll0 * (N - sum(dat[[4]])) + sum(llis)
    
    return(ans)
  }
  
  ## Vectorized version that exploits parallel computation
  llv <- function(l0, D, a, b) {
    future_pmap_dbl(
      .l = list(l0, D, a, b),
      .f = ll,
      .options = furrr_options(seed = TRUE)
    )
  }
}

## The log-likelihood for the experiment can be split in two elements the first
## for the cumulative density of the capture periods, the other one with the
## probability of death for those who are not trapped.

#' Log-likelihood for the capture periods (observed).
#'
#' @param loc Numeric 2-column matrix with trap coordinates.
#' @param pars List with parameters of the system: pds, diffusivity, alpha and
#'   beta.
#'
#' @return A function that takes the trap's id and the limits of a period and
#'  returns the log probability of be trapped in these conditions
#' @export
#'
#' @examples
logLik_obs <- function(loc, pars) {
  hi <- expected_specific_hazard_functions(loc, pars)
  s <- overall_survival(loc, pars, n_eval = 1e2)
  
  function(i, ts, te) {
    ni <- length(i)  # number of observations
    stopifnot(
      length(ts) == ni,
      length(te) == ni
    )
    
    ans <- numeric(ni)
    
    ## Iterate over each observation j
    for(j in seq.int(ni)) {
      f <- function(t) {
        hit <- hi(t)
        hit$h[hit$i == i[j]] * s(t)
      }
      ans[j] <- cubintegrate(
        f, ts[j], te[j],  relTol = 1e-2
      )$integral
    }
    
    return(log(ans))
  }
}

# logLik_non_obs(loc, pars)   # 6 - 10 min
#' Log-likelihood for the death case (unobserved).
#'
#' @param loc Numeric 2-column matrix with trap coordinates.
#' @param pars List with parameters of the system: pds, diffusivity, alpha and
#'   beta.
#'
#' @return A real value corresponding to the log probability of death in our
#' setting.
#' @export
#'
#' @examples
logLik_non_obs <- function(loc, pars) {
  
  s <- overall_survival(loc, pars, n_eval = 1e3)
  l0 <- -log(pars$pds)
  alter_s <- function(t){s(t/(1-t)) / (1 - t)**2}   # alternative version to integrate between 0 and 1
  int_s <- cubintegrate(
    alter_s, 0, 1, relTol = 1e-2
  )
  
  return(log(l0) + log(int_s$integral))
}

#' Simpson's rule integration
#'
#'
#' @param y vector of function evaluations in a sequence of an odd-number of
#'   evenly spaced points
#' @param h distance between evaluation points
#'
#' @return Integral of the function between the first and last points
#'
#' @examples
#' simpsons_rule(dnorm(seq(-1.96, 1.96, length = 11)), h = 2 * 1.96 / 10)
simpsons_rule <- function(y, h) {
  
  n <- length(y)
  
  ## n should be odd
  stopifnot(n %% 2 != 0, n > 2)
  
  coefs <- c(1, 4, rep(c(2, 4), times = (n - 3) / 2), 1)
  
  return(h * sum(coefs * y) / 3)
}

#' Sequence through a set of numbers
#'
#' @param through
#' @param length.out
#'
#' @return
#' @export
#'
#' @examples
#' length.out = 10; through = c(0, 1, 4, 8)
#' seq_interpolate(through, length.out)
seq_interpolate <- function(through, length.out) {
  through <- sort(through)
  rg = range(through)
  interval_lengths <- diff(through)
  lengths <- round(length.out * interval_lengths / sum(interval_lengths))
  
  ## Make sure lengths are even (for simpson's rule) and at least 4
  lengths <- pmax(4, lengths + lengths %% 2)
  
  ans <- pmap(
    list(
      from = head(through, n = -1),
      to = tail(through, n = -1),
      by = interval_lengths / lengths
    ),
    seq
  )
  
  # ans <- unique(do.call('c', sub_sequences))
  return(ans)
}



# VERIFICATIONS ---------------------------------------------------------------

## Numerical and graphical verifications about the simulation and the model.

#' Sequence of values around a point
#'
#' Return a numeric vector with a sequence of values around a central point.
#'
#' @param x Numeric. Central point
#' @param p Numeric. Fraction of the value subtracted and added to determine de
#'   limits of the sequence.
#' @param n Integer. Number of values to return. Use an odd number to find the
#'   original value in the sequence.
#'
#' @return Numeric vector.
#' @export
#'
#' @examples
#' seq_around(5, .1, 7)
#' seq_around(5, .2, 6)
#' seq_around(5, .3, 7)
seq_around <- function(x, p, n) {
  bounds <- x + c(-1, 1) * p * x
  seq(bounds[1], bounds[2], length.out = n)
}

#' Grid of values around a multidimensional point
#'
#' Return a data.frame with gridded values around a central point.
#'
#' @param x Numeric vector. Central point of any dimension. Possibly named.
#' @param p Numeric. Fraction of the value subtracted and added to determine de
#'   limits of the sequence.
#' @param n Integer. Number of values to return. Use an odd number to find the
#'   original value in the sequence.
#'
#' @return Data.frame
#' @export
#'
#' @examples
#' grid_around(c(D = 5, lambda = 2, alpha = .3), p = .1, n = 3)
grid_around <- function(x, p, n) {
  
  sequence_around <- function(x, p, n) {
    x * (1 + p * seq(-(n-1)/2, (n-1)/2))
  }
  
  x |>
    lapply(sequence_around, p, n) |>
    expand.grid() |>
    tibble()
  
}

#' Empirical density of individuals
#'
#' Retrieve a empirical density estimate of active individuals at given time
#' points from a simulation.
#'
#' @param obj A \code{sitModel} object.
#' @param t Positive number. Time (in days).
#' @param n Positive integer. Number of grid points in each direction. Can be
#'   scalar or a length-2 integer vector.
#' @param limits Numeric vector of length 4. (xmin, xmax, ymin, ymax). Calculated
#'   by default from the ranges of variation in the data.
#'
#' @return A tibble with variables t (time in days), x, y (coordinates) and z
#'   (empirical density estimate).
#' @export
#'
#' @examples
#' empirical_density(sit_scenario, t = c(5, 10, 15), n = 5)
empirical_density <- function(obj, t, n = 25, limits = NULL) {
  
  dens2df <- function(dens) {
    expand_grid(x = dens$x, y = dens$y) |>
      mutate(
        z = as.vector(dens$z)
      )
  }
  
  time_idx <- round(t / parms(obj)$tstep)
  
  loc <- out(obj)[time_idx] |>
    setNames(nm = t) |>
    bind_rows(.id = "t") |>
    tibble() |>
    select(t, x, y) |>
    mutate(t = as.numeric(t))
  
  if (is.null(limits)) limits <- c(range(loc$x), range(loc$y))
  
  ans <- loc |>
    nest(.by = t, .key = 'coord') |>
    mutate(
      emp_dens = map(
        coord,
        ~dens2df(
          MASS::kde2d(.$x, .$y, n = n, lims = limits)
        )
      )
    ) |>
    select(-coord) |>
    unnest(cols = c(t, emp_dens))
  
  # loc <- map(out(obj)[time_idx], tibble) |>
  #   map(~select(., x, y)) |>
  #   setNames(nm = t)
  #
  # loc <- map(out(obj)[time_idx], tibble) |>
  #   map(~select(., x, y)) |>
  #   setNames(nm = t)
  #
  # emp_dens <- map(loc, ~MASS::kde2d(.$x, .$y, n = n))
  #
  # ans <- map(emp_dens, dens2df) |>
  #   bind_rows(.id = "t") |>
  #   mutate(t = as.numeric(t))
  
  return(ans)
}

#' Density of individuals
#'
#' Plot the empirical spatial density of the active individuals at the given
#' time, and compare with the corresponding theoretical density.
#'
#' @param obj A \code{sitModel} object.
#' @param t Positive number. Time (in days).
#' @param n Positive integer. Number of grid points in each direction. Can be
#'   scalar or a length-2 integer vector.
#'
#' @export
#'
#' @examples
#' empirical_vs_theoretical_density(sit_scenario_dispersion_only, t = 10)
#' empirical_vs_theoretical_density(sit_scenario, t = 10)
empirical_vs_theoretical_density <- function(obj, t, n = 25) {
  
  
  var_theo <- function(t) 2 * parms(obj)$diffusivity * t
  
  dat <- empirical_density(obj, t, n) |>
    rename(z_Empirical = z) |>
    mutate(
      z_Theoretical = dnorm(x, sd = sqrt(var_theo(t))) * dnorm(y, sd = sqrt(var_theo(t)))
    ) |>
    pivot_longer(
      cols = starts_with("z_"),
      names_prefix = "z_",
      names_to = "type",
      values_to = "z"
    ) |>
    mutate(t = paste("Day", t))
  
  ## Clean up environment shit that will be stored in ggplot layers.
  rm(obj)
  
  ans <- ggplot(dat) +
    aes(x, y) +
    facet_grid(t ~ type) +
    geom_raster(
      aes(fill = z),
      show.legend = FALSE
    ) +
    scale_fill_viridis_c() +
    coord_fixed() +
    labs(
      title = paste("Density of individuals"),
      x = NULL,
      y = NULL
    ) +
    theme_ipsum(grid = "") +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    )
  
  ## ans contains a copy of the environment, which includes the entire
  ## obj object, which is huge. Remove this prior to returning the result.
  ## https://stackoverflow.com/a/77373906
  ans$plot_env <- rlang::new_environment()
  
  return(ans)
}




#' Estimates of variance
#'
#' Plot scatter plot of the empirical vs theoretical estimates of variance
#' of the individuals at each time step of the simulation.
#'
#' @param obj A \code{sitModel} object.
#'
#' @export
#'
#' @examples
#' estimates_variance(sit_scenario_dispersion_only)
estimates_variance <- function(obj) {
  
  dat <- dispersion_behaviour(obj) |>
    select(the_var, emp_var)
  
  ## Clean up environment shit that will be stored in ggplot layers.
  rm(obj)
  
  ans <- ggplot(dat) +
    aes(the_var, emp_var) +
    geom_abline(intercept = 0, slope = 1, colour = 'grey') +
    geom_line() +
    labs(
      x = "Theoretical variance",
      y = "Empirical variance",
    ) +
    theme_ipsum(grid = "")
  
  ## ans contains a copy of the environment, which includes the entire
  ## obj object, which is huge. Remove this prior to returning the result.
  ## https://stackoverflow.com/a/77373906
  ans$plot_env <- rlang::new_environment()
  
  return(ans)
}


#' Estimates of diffusivity
#'
#' Plot an histogram of the empirical estimates of diffusivity from the
#' empirical dispersion of the individuals at each time step of the simulation.
#'
#' @param obj A \code{sitModel} object.
#'
#' @export
#'
#' @examples
#' estimates_diffusivity(sit_scenario_dispersion_only)
estimates_diffusivity <- function(obj) {
  
  dat <- dispersion_behaviour(obj) |>
    filter(est_D > 0) |>
    select(est_D)
  
  xi <- parms(obj)$diffusivity
  
  ## Clean up environment shit that will be stored in ggplot layers.
  rm(obj)
  
  ans <- ggplot(dat) +
    aes(est_D) +
    geom_histogram(bins = 21) +
    geom_vline(
      xintercept = xi,
      colour = 'darkorange'
    ) +
    labs(x = "Diffusivity", y = NULL) +
    theme_ipsum(grid = "") +
    theme(
      axis.text.y = element_blank()
    )
  
  ## ans contains a copy of the environment, which includes the entire
  ## obj object, which is huge. Remove this prior to returning the result.
  ## https://stackoverflow.com/a/77373906
  ans$plot_env <- rlang::new_environment()
  
  return(ans)
}


#' Dispersion behaviour
#'
#' Retrieve summaries of the dispersion over time from a simulation object.
#'
#' @param obj A \code{sitModel} object.
#'
#' The empirical estimate of diffusivity is isolated from Var = 2tD, using
#' the empirical variance.
#'
#' @return A tibble with the empirical and theoretical variance of the dispersion
#' and an the corresponding empirical estimate of the diffusivity.
#'
#' @export
#'
#' @examples
dispersion_behaviour <- function(obj) {
  
  out(obj) |>
    enframe(name = 'id', value = 'position') |>
    mutate(
      t = do.call(seq, as.list(times(obj))),
      coord = map(position, ~ select(., x, y) |> as.matrix()),
      emp_var = map_dbl(coord, ~var(.) |> diag() |> mean()),
      the_var = 2 * t * parms(obj)$diffusivity,
      est_D = emp_var / t / 2
    )
  
}


#' Empirical survival function
#'
#' Fraction of individuals that are still active (alive and not caught) as
#' a function of time.
#'
#' @param obj A \code{sitModel} object.
#'
#' @export
#'
#' @examples
#' survival(sit_scenario)
survival <- function(obj, plot = TRUE) {
  
  dat <- sim_results(obj) |>
    count(t) |>
    mutate(p = n / n[1])
  
  attr(dat, "trap_counts") <- NULL
  
  if (!plot) return(dat)
  
  ans <- dat |>
    ggplot(aes(t, p)) +
    geom_line() +
    labs(
      x = "Time (days)",
      y = "Empirical survival function"
    ) +
    theme_ipsum(grid = "Y")
  
  ## ans contains a copy of the environment, which includes the entire
  ## obj object, which is huge. Remove this prior to returning the result.
  ## https://stackoverflow.com/a/77373906
  ans$plot_env <- rlang::new_environment()
  
  return(ans)
}


#' Extract summaries of dispersion behaviour from replicated diffusion process
dispersion_replicates <- function(x) {
  map(x, dispersion_behaviour) |>
    bind_rows(.id = "rep") |>
    select(-coord)  # Halves the object's size
}


#' Plot a simulation at a time step
#'
#' @param obj Simulation results. Output from `sim_results()`.
#' @param time Positive number. Target time frame to export.
#' @param traps Data frame with trap locations and counts. Attribute
#'   'trap_counts' from the output of `sim_results`.
#' @param fn Character. File name for the exported figure.
#'
#' @export
#'
#' @examples
#' res <- sim_results(sit_scenario)
#' plot_simulation_step(
#'    res,
#'    time = 15,
#'    traps = attr(res, 'trap_counts')
#' )
plot_simulation_step <- function(obj, time, traps) {
  
  time_steps <- tibble(
    day = time
  )
  
  inner_join(
    time_steps,
    obj,
    by = join_by(day == t)
  ) |>
    ggplot(aes(x, y)) +
    facet_wrap(~ day, labeller = "label_both") +
    geom_point(
      size = .4,
      alpha = 0.7, show.legend = FALSE, na.rm = TRUE
    ) +
    geom_point(
      data = inner_join(
        time_steps,
        traps |>
          left_join(
            attr(traps, 'loc'),
            by = 'id'
          ),
        by = join_by(day == t)
      ),
      aes(size = n_cum),
      col = 'darkorange',
      shape = 15
    ) +
    coord_fixed() +
    # theme_void() +
    labs(
      size = "Catch"
    ) +
    theme_void()
  
}


#' Plot log-likelihood
plot_loglik <- function(x, lower_cutoff = -Inf) {
  
  x |>
    ## Put a lower bound to the likelihood to avoid obscuring the
    ## interesting part around the maximum
    mutate(ll = pmax(ll, lower_cutoff)) |>
    ggplot() +
    facet_nested(
      "Capture rate (a)" + round(a, 3) ~ "Capture extent (b)" + round(b, 3),
      # labeller = "label_both"
    ) +
    aes(l0, D, fill = ll) +
    geom_raster() +
    geom_point(
      data = ~ filter(., ll >= max(ll, na.rm = TRUE) - 2),
      size = 2,
      colour = "blue"
    ) +
    labs(
      x = "Death rate (l0)",
      y = "Diffusivity (D)",
      fill = "Log-likelihood"
    ) +
    scale_x_log10() +
    scale_y_log10() +
    scale_fill_viridis_c(option = "B")
}


#' Give a cue to run a target on a specific list of nodes
#'
#' @examples
#' # use in the argument cue of tar_target() to ensure that the target
#' # is only run in the computer named "fmcd"
#' run_in_node("fmcd")
run_in_node <- function(...) {
  target_nodes <- list(...)
  current_node <- Sys.info()[["nodename"]]
  do_run <- isTRUE(current_node %in% target_nodes)
  ifelse(do_run, "always", "never")
}



#' Simulation of a Brownian Motion
#'
#' Simulate the trajectory of a particle following a Brownian motion.
#' This is, at each time step, the displacement in x and y are simulated
#' from a Gaussian distribution.
#'
#' We use a standard Gaussian for this simulation. An increased diffusivity
#' can be obtained by increasing the number of steps in the same period of time.
#'
#' @param n_steps Number of time steps in the simulation.
#' @param n_days Number of days in total (i.e. fixing the duration of the time step).
#' @param seed Arbitrary number for starting the random number generator.
#'
#' @return
#' @export
#'
#' @examples
#' brownian_motion(10, 2, 1)
brownian_motion <- function(n_steps, n_days, seed) {
  
  set.seed(seed)
  data.frame(
    t = seq.int(n_steps) / n_steps * n_days,
    x = cumsum(rnorm(n_steps)),
    y = cumsum(rnorm(n_steps))
  ) |>
    mutate(
      day = cut(t, seq.int(n_days + 1) - 1, labels = seq.int(n_days))
    )
  
}

#' Computation of approximate and empirical of capture probability
#'
#' @param obj Scenario created with simulate_point_release()
#'
#' @return A tibble with the theoretical and empirical capture probability
#' with the corresponding confidence interval (95%) and the quantile of the
#' model.
#' @export
#'
#' @examples
#' capture_prob_data(sit_scenario)
capture_prob_data <- function(obj){
  
  N <- nrow(init(obj))
  
  # Empirical trap counts
  captures <- trap_counts(obj)
  total_captures <- sum(captures$n)
  
  # Simulations Parameters
  pars <- simulation_parameters(obj)
  loc <- trap_coords(obj) |> select(-id) |> as.matrix()
  
  # Survival end specific hazard functions
  S <- overall_survival(loc, pars)
  h <- expected_specific_hazard_functions(loc, pars)
  
  ## Partition of the simulation time span into 200 time intervals
  max_time <- ceiling(max(capture_times(obj)$t))
  time_step <- max_time / 200
  time_points <- seq(0, max_time, time_step)
  
  
  # Theoretical probability approximation
  # Approximate probability of termination at each time-step of 0.1, for each
  # cause.
  dat <- h(time_points) |>
    rowwise() |>
    mutate(z = h * S(t))
  
  summary_dat <-
    dat |>
    ## Probability of termination cause:
    ## Sum of probabilities of termination at each time-step during the full
    ## span of the experiment.
    group_by(i) |>
    summarise(
      prob_theoretical = time_step * sum(z),
    ) |>
    ## Empirical probability estimation
    left_join(
      captures |>
        group_by(id) |>
        summarise(
          sum_n = sum(n),
          freq = sum_n / N
        ) |>
        add_row(id = 0L, freq = 1 - total_captures / N) |>
        mutate(i = factor(id)) |>
        select(-id),
      by = join_by('i')
    ) |>
    ## Confidence intervals
    mutate(
      inf = pmax(0, freq - 1.96 * sqrt(freq*(1-freq) / N)),
      sup = pmax(0, freq + 1.96 * sqrt(freq*(1-freq) / N)),
    ) |>
    ## Quantile for the theoretical model
    mutate(
      quant_inf = qbinom(0.025, N, prob_theoretical),
      quant_sup = qbinom(0.975, N, prob_theoretical)
    )
  
  return(summary_dat)
  
}

# OTHERS ---------------------------------------------------------------------

## Documented functions which exist but not necesserarly used in the target file.

#' Distribution parameters of termination causes under fixed position.
#'
#' Return the probabilities of termination causes, where the first position
#' corresponds to death and the remaining to the traps, in order.
#'
#' @param x Numeric vector. Distance to traps.
#' @param pds Number between 0 and 1. Probability of daily survival.
#' @param a Positive number. Daily capture rate at distance 0.
#' @param b Positive number. Decay rate with respect to squared distance.
#'
#' @return Numeric vector, one element longer than x. Probabilities of death
#' and capture in each trap.
#'
#' @examples
#' probs_causes_fixed(x = c(1, 2, 4), pds = .8, a = .8, b = .2)
probs_causes_fixed <- function(x, pds, a, b) {
  
  l0 <- -log(pds)
  li <- exp(log(a) - b * x ** 2)
  
  lambdas <- c(l0, li)
  L <- sum(lambdas)
  
  return(lambdas / L)
}

hazard_rates <- function(x) {
  
  ## Position and directions of the active individuals by time
  res <- sim_results(x)
  
  ## Trap counts by time
  trap_counts <- attr(res, 'trap_counts') |>
    left_join(
      parms(x)$traps.loc,
      by = 'id'
    )
  
  ## Total released individuals
  N <- parms(x)$ninds
  
  ## Extrapolate the total number of catches by trap
  ## to account for the remaining active individuals
  specific_totals <- trap_counts |>
    group_by(x, y) |>
    filter(t == max(t)) |>
    ungroup() |>
    mutate(
      ## Expected number of catches by trap, if the proportion is maintained
      N = n_cum / sum(n_cum) * N
    ) |>
    select(id, x, y, N)
  
  specific_h <- trap_counts |>
    left_join(
      specific_totals,
      by = join_by(id, x, y)
    ) |>
    mutate(
      f = n / N,
      s = (N - n_cum) / N,
      h = f / s / parms(.env$x)$tstep  # Scale to daily rate
    )
  
  overall_h <- trap_counts |>
    group_by(t) |>
    summarise(
      n = sum(n),
      n_cum = sum(n_cum)
    ) |>
    mutate(
      id = 0,
      f = n / N,
      s = (N - n_cum) / N,
      h = f / s / parms(.env$x)$tstep  # Scale to daily rate
    )
  
  ans <- bind_rows(
    specific_h |>
      ungroup() |>
      select(id, t, n, n_cum, f, s, h),
    overall_h |>
      select(id, t, n, n_cum, f, s, h)
  ) |>
    filter(t > 0) |>
    left_join(
      ## Add column with distance of trap to origin
      parms(x)$traps.loc |>
        mutate(
          z = sqrt(x**2 + y**2)
        ) |>
        select(id, z),
      by = 'id'
    )
  
  return(ans)
}


#' Error function
#'
#' The error or Phi function is a variant of the cumulative normal (or Gaussian)
#' distribution.
#'
#' This implementation is borrowed from the package {pracma}.
#'
#' @param x Vector of real numbers
#'
#' @return Real numbers, between -1 and 1.
#' @export
#'
#' @examples
#' erf(1); 2*pnorm(sqrt(2)*1) -1
erf <- function(x) {
  pchisq(2 * x^2, 1) * sign(x)
}

#' Cumulative probability of capture
#'
#' Given parameters a, b and c, characterising the location and time of the
#' survey, return the cumulative probability of capture from release to survey
#' time.
cpit <- function(x, a, b, c) {
  a / 2 *
    sqrt(pi / b) *
    exp(-2 * sqrt(b * c)) *
    (
      -erf((sqrt(c) - sqrt(b) * x) / sqrt(x)) +
        exp(4 * sqrt(b * c)) *
        (erf((sqrt(c) + sqrt(b) * x) / sqrt(x)) - 1) + 1
    )
}
#' Daily probability of capture
#'
#' Given parameters a, b and c, characterising the location and time of the
#' survey, return the daily probability of capture at day x
#' @examples
#' pit(8:14, a = .00163, b = 0.223, c = 28.9)
pit <- function(x, a, b, c) {
  ans <- a * x ^ (-1 / 2) * exp(-b * x - c / x)
  ans[is.nan(ans)] <- 0
  return(ans)
}

#' Integer Check
#'
#' Checks whether the submitted values are whole numbers.
#'
#' @param x A `numeric` vector
#'
#' @return
#' A `logical` vector of the same size indicating whether the value is an
#' `whole` number or not.
#'
#' @export
#' @examples
#' is_whole(2.3)
#' is_whole(-3)
#' is_whole(4)
#' is_whole(c(1,2,3))
#' is_whole(c(.4,.5,.6))
#' is_whole(c(7,.8,9))
is_whole = function(x, tol = sqrt(.Machine$double.eps)) {
  if (!is.numeric(x)) return(FALSE)
  vapply(x - round(x), \(.) isTRUE(all.equal(., 0, tolerance = tol)), NA)
}

#' Specific hazard rate
#'
#' This is the hazard rate in the presence of a single termination cause
#' (e.g. capture in a trap, without mortality or other traps).
#'
#' The individual is assumed to be released in the origin of coordinates at time
#' 0.
#'
#' @param t Positive number. Time where the hazard is evaluated.
#' @param D Positive number. Diffusion parameter.
#' @param beta Positive number. Decay parameter of the hazard.
#' @param alpha Positive number. Hazard at distance 0.
#' @param x0 Numeric vector. Coordinates of the trap.
#'
#' @return A positive number. Hazard rate of an individual at time t.
#'
#' @example
#' hazard_one(1, x0=c(1,1))
#' hazard_one(1:10, x0=c(1,1))
#'
hazard_one <- function(t, D=2, alpha=0.8, beta=0.2, x0){
  
  n <- sum(x0**2)  # euclidean norm of x0 squared
  
  c <- 1 / (1 + 4 * D * t * beta)
  
  return(alpha * exp(-beta * n * c) * c)
}


#' Capture pdf in a single-trap setting.
#'
#' @param t Positive number. Time where the hazard is evaluated.
#' @param D Positive number. Diffusion parameter.
#' @param beta Positive number. Decay parameter of the hazard.
#' @param alpha Positive number. Hazard at distance 0.
#' @param x0 Numeric vector. Coordinates of the trap.
#'
#' @returnA positive number. Pdf at time t.
#' @export
#'
#' @examples
#' density_one(0.5,2,.8,.2,c(1,1))
#' density_one(1:10,2,.8,.2,c(1,1))
density_one <- function(t, D=2, alpha=0.8, beta=0.2, x0){
  
  h_uni <- function(t, D, alpha, beta, x0){
    integrate(\(.) hazard_one(.,D,alpha,beta,x0),0,t)$value
  }
  
  h <- Vectorize(h_uni, vectorize.args = 't')
  
  return(hazard_one(t, D, alpha, beta, x0) * exp(-h(t, D, alpha, beta, x0)))
}



#' Cumulative hazard function using an alternative form of the integral.
#'
#' @param t Positive number. Time where the hazard is evaluated.
#' @param D Positive number. Diffusion parameter.
#' @param beta Positive number. Decay parameter of the hazard.
#' @param alpha Positive number. Hazard at distance 0.
#' @param x0 Numeric vector. Coordinates of the trap.
#'
#' @return A list of two elements :
#'      - time : which is referred to the vector t.
#'      - cumhaz : which is the values of the cumulative hazard function at
#'      time t.
#'
#' @export
#'
#' @examples
cumhazard_one <- function(t,  D=2, alpha=0.8, beta=0.2, x0){
  # Computation of the norm
  r <- as.numeric(t(x0) %*% x0)
  
  # Definition of the alternative integral (see calculations.html)
  h <- function(u){
    (1 / u) * exp(-beta * r * u)
  }
  
  tt <- length(t)
  
  H <- rep(0, tt)
  
  H[1] <- integrate(h, lower = 1 / (1 + 4 * D * beta * t[1]),
                    upper = 1, abs.tol = 1e-4)$value
  
  for(i in 2:tt){
    H[i] <- H[i - 1] + integrate(h, lower = 1 / (1 + 4 * D * beta * t[i]),
                                 upper = 1 / (1 + 4 * D * beta * t[i - 1]),
                                 abs.tol = 1e-4)$value
  }
  
  return(list(time = t, cumhaz = (alpha / (4 * D * beta)) * H))
}








#' Lambdas at a position x, with multiple traps.
#'
#' @param x Vector of size 2 with the position of the individual.
#' @param A Trap coordinates.
#' @param alpha Parameter hazard rate.
#' @param beta Parameter hazard rate.
#'
#' @return Numeric vector of size equal to the number of traps.
#' @export
#'
#' @examples
#' A <- matrix(
#'   c(5, 10, 20, 0, 5, 5, 0, 40),
#'   ncol = 2
#' )
#' lambda_x(c(1, 2), A, .8, .2)
#' lambda_x(c(1, 2), A[1, ], .8, .2)
lambda_x <- function(x, A, alpha=0.8, beta=0.2){
  
  if (!is.matrix(A)) A <- t(A)
  
  A_shift <- t(t(A) - x)
  
  r_x <- as.numeric(A_shift ** 2 %*% c(1, 1))
  
  return(alpha * exp(-beta * r_x))
}

#' Monte Carlo computation of the expected value of the probability of capture
#' in a trap
#'
#' @param t Numeric vector of time points.
#' @param D Diffusion.
#' @param A Locations traps.
#' @param i Trap number.
#' @param alpha Parameter dispersion.
#' @param beta Parameter dispersion.
#' @param n_iter Monte Carlo iterations.
#'
#' @examples
#' A <- matrix(
#'   c(5, 10, 20, 0, 5, 5, 0, 40),
#'   ncol= 2
#' )
#' exp_mc(5:10, D = 2, A, i = 2, alpha = .8, beta = .2, n_iter = 100)
#'
exp_mc_single <- function(t, D=2, A, i=1, alpha=.8,
                          beta=.2, n_iter=1000){
  
  val <- 0
  
  gamma_t <- (1 + 8 * D * beta * t) / (8 * D * t)
  
  Y <- cbind(
    rnorm(n_iter, beta*A[i,1]/gamma_t, 1 / (2 * gamma_t)),
    rnorm(n_iter, beta*A[i,2]/gamma_t, 1 / (2 * gamma_t))
  )
  
  for(k in 1:n_iter) {
    l_x <- lambda_x(Y[k,], A, alpha, beta)
    val <- val + l_x[i]/sum(l_x)
  }
  return(val/n_iter)
}

exp_mc <- Vectorize(exp_mc_single,
                    vectorize.args = 't')

#' Capture hazard in a multiple-trap setting
#'
#'
#' @param t Positive number. Time where the hazard is evaluated.
#' @param D Positive number. Diffusion parameter.
#' @param beta Positive number. Decay parameter of the hazard.
#' @param alpha Positive number. Hazard at distance 0.
#' @param k Positive integer. Trap number.
#' @param A 2-column matrix with coordinates of traps with respect to the
#'   release point, in some spatial units (e.g. m).
#' @param mc.iter Positive integer. Number of Monte Carlo iterations for the
#'   integration step.
#'
#' @return A positive number. Hazard rate of an individual at time t for trap k.
#'
#' @example
#' A <- matrix(
#'   c(5, 10, 20, 0, 5, 5, 0, 40),
#'   ncol = 2
#' )
#' hazard_multi(5, D = 2, 0.8, 0.2, A, k = 2, mc.iter = 200)
#' hazard_multi(1:10, D = 2, 0.8, 0.2, A, k = 2, mc.iter = 200)
hazard_multi <- function(t, D, alpha, beta, A, k, mc.iter = 2000) {
  
  a <- A[k, ]
  
  ans <- hazard_one(t, D, alpha, beta, a) * exp_mc(t, D, A, k, alpha, beta, mc.iter)
  return(as.numeric(ans))
}


#' Capture pdf in a multi-trap setting.
#'
#' @param t Numeric vector. Times where the pdf is evaluated.
#' @param D Diffusion parameter.
#' @param D Positive number. Diffusion parameter.
#' @param beta Positive number. Decay parameter of the hazard.
#' @param alpha Positive number. Hazard at distance 0.
#' @param k Positive integer. Trap number.
#' @param A 2-column matrix with coordinates of traps with respect to the
#'   release point, in some spatial units (e.g. m).
#' @param mc.iter Positive integer. Number of Monte Carlo iterations for the
#'   integration step.
#'
#' @return  A numeric vector. Pdf at time t for trap k.
#' @export
#'
#' @examples
#' A <- matrix(
#'   c(5, 10, 20, 0, 5, 5, 0, 40),
#'   ncol = 2
#' )
#' density_multi(1:10, D = 2, 0.8, 0.2, A, k = 1, mc.iter = 200)
#'
density_multi <- function(t, D=2, alpha=0.8, beta=0.2, A, k, mc.iter = 2000){
  
  f <- \(.) hazard_multi(., D, alpha, beta, A, k, mc.iter)
  
  n <- length(t)
  
  h <- rep(0, n)
  
  for(i in 2:n){
    h[i] <- h[i-1] + integrate(f, t[i-1], t[i], rel.tol = 1e-3)$value
  }
  return(hazard_multi(t, D, alpha, beta, A, k, mc.iter) * exp(-h))
}

# Version avec approximation de l'espérance pour un t assez proche (bien plus rapide pour les grandes valeurs)
density_multi2 <- function(t, D=2, alpha=0.8, beta=0.2, A, k, mc.iter = 2000, tol = 1e-2){
  
  n <- length(t)
  
  error <- abs(1 - (8 * D * beta * t) / (1 + 8 * D * beta * t))
  
  h <- rep(0, n)
  
  # Approached hazard function with approximation of the true expectation
  f <- \(.) hazard_multi(., D, alpha, beta, A, k, mc.iter)
  
  ## Approached hazard function with alternative expectation
  exp_approx <- exp_mc(max(t), D, A, k, alpha, beta, mc.iter)
  
  a <- A[k, ]
  
  f_approx <- \(.) exp_approx*hazard_one(., D, alpha, beta, a)
  
  for(i in 2:n){
    # if the error is low enough, we use the alternative hazard
    if(error[i]>tol){
      h[i] <- h[i-1] + integrate(f, t[i-1], t[i], rel.tol = 1e-3)$value
    }else{
      h[i] <- h[i-1] + integrate(f_approx, t[i-1], t[i], rel.tol = 1e-3)$value
    }
    
  }
  return(hazard_multi(t, D, alpha, beta, A, k, mc.iter) * exp(-h))
}


#' Capture times and traps conditional on capture.
#'
#' Simulates a random variable from a hyperexponential distribution.
#'
#' @param n Number of values to simulate.
#' @param lambda Vector of positive numbers. Exponential parameters.
#' @param p Mixture probabilities.
#'
#' @return A list with capture times and capture trap.
#' @examples
#' rhypexp(10, c(2, 3, 7), p = c(.3, .5, .2))
rhypexp <- function(n, lambda, p){
  
  I <- length(p)
  
  traps <- sample(1:I, n, replace = TRUE, prob = p)
  
  l <- lambda[traps]
  
  res <- rep(0, n)
  
  for(i in 1:n){
    res[i] <- rexp(1, l[i])
  }
  return(list("simulation" = res, "traps" = traps))
}



#' Density of capture time in presence of multiple traps.
#'
#' Capture time in no matter which trap, conditional on capture.
#'
#' @param x Numeric vector. Times of evaluation of density.
#' @param lambda Vector of positive numbers. Exponential parameters.
#' @param p Mixture probabilities.
#'
#' @examples
#' dhyperexp(0.3, lambda = c(2, 3, 7), p = c(.3, .5, .2))
#' dhyperexp(x = seq(0, 10, by = 0.1), lambda = c(2, 3, 17), p = c(.3, .5, .2))
dhyperexp <- function(x, lambda, p) {
  
  lambdax <- outer(lambda, x)
  
  ans <- p %*% (diag(lambda) %*% exp(-lambdax))
  
  return(as.numeric(ans))
}

#' Global hazard rate of capture
#'
#' Risk of capture in any of multiple traps as a function of time.
#'
#' @param x Numeric vector. Times of evaluation of density.
#' @param lambda Vector of positive numbers. Exponential parameters for each trap.
#' @param p Mixture probabilities.
#'
#' @examples
#' hhyperexp(0.3, lambda = c(2, 3, 7), p = c(.3, .5, .2))
#' hhyperexp(x = seq(0, 10, by = 0.1), lambda = c(2, 3, 17), p = c(.3, .5, .2))
hhyperexp <- function(x, lambda, p) {
  
  lambdax <- outer(lambda, x)
  
  denom <- p %*% exp(-lambdax)
  
  ans <- dhyperexp(x, lambda, p) / denom
  
  return(as.numeric(ans))
}

#' Simulation with Brownian motion
#'
#' @param h Time step.
#' @param m_max Maximum time of simulation.
#' @param D Diffusion parameter.
#'
#' @examples
#' brownian(1, 15, 2)
#'
brownian <- function(h, m_max, D){
  
  pas1 <- rnorm(m_max, 0, sqrt(4 * D * h))
  
  pas2 <- rnorm(m_max,0,sqrt(4 * D * h))
  
  return(cbind(cumsum(pas1), cumsum(pas2)))
}


#' Alternative computation of the cumulative hazard, with approximation.
#'
#' @param t Positive value : time.
#' @param D Positive parameter : diffusion
#' @param alpha Positive number. Linked with the capture risk at 0 meter.
#' @param beta Positive number. Linked with alpha and the capture risk at 5 meters.
#' @param x0 Numeric vector : trap location.
#'
#' @return Compute the value of the cumulative hazard function with an approximation
#' in the neighborhood of infinity.
#' @export
#'
#' @examples
#'
#'
H <- function(t,  D=2, alpha=0.8, beta=0.2, x0){
  # Computation of the norm
  r <- as.numeric(t(x0) %*% x0)
  
  # Definition of the alternative integral (see calculations.html)
  h <- function(u){
    (1 / u) * exp(-beta * r * u)
  }
  
  H <- 0
  
  # Bound calculation for the log-approximation
  delta <- -log(0.9995)
  
  t1 <- (beta * r - delta) / (delta * 4 * D * beta)
  
  if(t<t1){
    H <- cubintegrate(h, lower = 1 / (1 + 4 * D * beta * t),
                      upper = 1, method = "pcubature")$integral
  }else{
    H <-
      cubintegrate(h, lower = 1 / (1 + 4 * D * beta * t1), upper = 1,
                   method = "pcubature")$integral
    +
      (log(1+4*D*beta*t) - log(1+4*D*beta*t1))
  }
  
  return(alpha / (4 * D * beta) * H)
}

# Vectorized version for the t variable
cumhazard_one2 <- Vectorize(H, vectorize.args = 't')

#' Density without the normalization constant.
#'
#' @param t Positive value : time.
#' @param xk Numeric vector : trap location.
#' @param D Positive parameter : diffusion.
#' @param alpha Positive number. Linked with the capture risk at 0 meter.
#' @param beta Positive number. Linked with alpha and the capture risk at 5 meters.
#' @param lambda Positive number : survival parameter.
#'
#' @return Compute the density of observed times without the normalization constant.
#' @export
#'
#' @examples
#'
#'
s <- function(t, xk, D, alpha, beta, lambda){
  # Density of capture
  H1 <- cumhazard_one2(t,
                       D = D,
                       alpha = alpha,
                       x0 = xk,
                       beta = beta)
  
  h1 <- hazard_one(t,
                   D = D,
                   alpha = alpha,
                   x0 = xk,
                   beta = beta)
  
  # Survival probability
  s1 <- (1 - pexp(t, lambda))
  
  return(s1 * h1 * exp(-H1))
}


#' Specific log-likelihood for one trap.
#'
#' @param data Tibble. Output from the simulation.
#' @param x Numeric vector : trap position
#' @param D Positive parameter : diffusion.
#' @param alpha Positive number. Linked with the capture risk at 0 meter.
#' @param beta Positive number. Linked with alpha and the capture risk at 5 meters.
#' @param lambda Positive number : survival parameter.
#'
#' @return Compute the value of the log-likelihood for one trap setting.
#' @export
#'
#' @examples
#'
#'
pllh <- function(data, x, D, alpha, beta, lambda){
  # Capture time building
  capture <- inverse.rle(
    structure(
      list(
        lengths = data$n,
        values = (data$t - 1) * params1$tstep
      ),
      class = 'rle'
    )
  )
  
  sj  <- \(.) s(., x, D, alpha, beta, lambda)
  
  result <- sum(log(sj(capture) / integrate(sj, 0, Inf)$value))
  
  return(result)
}

## Proof of concept:
##
## Target factory to call a function that makes a plot and saves it
## in a associated file.
##
## See documentation at
## https://wlandau.github.io/targetopia/contributing.html#target-factories
##
## I get an error when calling the function that makes the plot
## the created factory does not find the function.
tar_save_figure <- function(fun, width, height, dir = "reports") {
  
  fun_name = deparse1(substitute(fun))
  
  ## remove calling parentheses for naming objects
  fun_label <- gsub("\\(.*\\)$", "", fun_name)
  
  file_name = paste(gsub("^plot_", "", fun_label), "png", sep = ".")
  
  target_name = gsub("^plot_", "fig_", fun_label)
  
  fn <- file.path(dir, file_name)
  
  command <- substitute(
    {
      ggsave(
        filename = fn,
        plot = fun,
        width = width,
        height = height,
        bg = "white"
      )
      return(fn)
    },
    env = list(
      fn = fn,
      fun = as.symbol(fun_name),
      width = width,
      height = height
    )
  )
  list(
    tar_target_raw(
      target_name,
      command,
      format = "file"
    )
  )
  
}

