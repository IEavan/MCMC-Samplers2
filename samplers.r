ID <-  c(1,6,3,0,0,2,2)
density <- function(x) {
    t1 <- ID[7] * exp(-sin((ID[1] * x ** 2)/(15 - ID[1])) - (x - 3 - ID[2] * pi) ** 2 / (2 * (5 + ID[3]) ** 2))
    t2 <- 2 * (1 + ID[7]) * exp(- x**2 /32)
    t3 <- (10 - ID[7]) * exp(-cos(ID[4] * x ** 2 / (15 + ID[4])) - (x + 3 + ID[5] * pi)**2/(2*(5+ID[6])**2))
    return(t1 + t2 + t3)
}

gmm_bound <- function(x) {
    t1 <- ID[7] * exp(1 - (x - 3 - ID[2] * pi) ** 2 / (2 * (5 + ID[3]) ** 2))
    t2 <- 2 * (1 + ID[7]) * exp(- x**2 /32)
    t3 <- (10 - ID[7]) * exp(-1 - (x + 3 + ID[5] * pi)**2/(2*(5+ID[6])**2))
    return(t1 + t2 + t3)
}

MH <- function(num_samples, x0 = 0, temperature_scheme = rep(1, num_samples), sd = 1){
    samples <- matrix(nrow = num_samples + 1)
    samples[1] <- x0
    for (i in 1:num_samples) {
        proposal <- rnorm(1, mean = samples[i], sd = sd)
        if (runif(1) <= (density(proposal) / density(samples[i])) ** temperature_scheme[i]) {
            samples[i+1] <- proposal
        } else {samples[i+1] <- samples[i]}
    }
    return(samples[2:num_samples+1])
}

simulated_annealing <- function(num_samples, max_temp = 10, x0 = 0, sd = 1) {
    temperature_scheme = seq(1, max_temp, length.out = num_samples)
    samples <- MH(num_samples, x0 = x0, temperature_scheme = temperature_scheme, sd = sd)
    
    mode <- NULL
    density_mode <- -1
    for (sample in samples){
        if (density(sample) > density_mode){
            mode <- sample
            density_mode <- density(sample)
        }
    }
    return(mode)
}

slice_sampler <- function(num_samples, x0 = 0, y0 = 1) {
    samples <- matrix(nrow = num_samples + 1, ncol = 2)
    samples[1,] <- c(x0,y0)
    for (i in 1:num_samples){
        lskew <- sqrt(-130 * (log(samples[i,2]) - 10))
        rskew <- sqrt(-128 * (log(samples[i,2]) - 2 * exp(1)))
        min_x = min(0.35 - lskew, 3 + 6 * pi - rskew)
        max_x = max(0.35 + lskew, 3 + 6 * pi + rskew)
        
        new_x <- NULL
        while (TRUE) {
            proposal <- runif(1, min_x, max_x)
            if (density(proposal) > samples[i,2]) {
                new_x <- proposal
                break
            }
        }
        
        new_y <- runif(1, 0, density(new_x))
        samples[i+1,] <- c(new_x, new_y)
    }
    return(samples[2:num_samples+1,])
}

rgmm <- function(num = 1) {
    choice <- runif(num)
    samples <- c(rnorm(length(choice[choice < 0.4937]), mean = 3 + 6 * pi, sd = 8))
    samples <- c(samples, rnorm(length(choice[0.4937 <= choice & choice < 0.7661]), mean = 0, sd = 4))
    samples <- c(samples, rnorm(length(choice[0.7661 <= choice]), mean = -3, sd = 7))
    return(samples)
}

dgmm <- function(x) {
    return(0.4937 * dnorm(x, mean = 3 + 6 * pi, sd = 8) +
           0.2724 * dnorm(x, mean = 0, sd = 4) +
           0.2339 * dnorm(x, mean = -3, sd = 7))
}

indep_sampler <- function(num_samples, x0 = 0){
    samples <- matrix(nrow = num_samples + 1)
    samples[1] <- x0
    for (i in 1:num_samples) {
        proposal <- rgmm(1)
        acceptance_prob <- (density(proposal) / density(samples[i])) * (dgmm(samples[i]) / dgmm(proposal))
        if (runif(1) < acceptance_prob) {
            samples[i + 1] <- proposal
        } else {
            samples[i + 1] <- samples[i]
        }
    }
    return(samples[2:num_samples+1])
}

rej_sampler <- function(num_samples = 1){
    samples <- c()
    pseudo_m <- 1
    while (length(samples) < num_samples) {
        quantity <- pseudo_m * (num_samples - length(samples))
        proposals <- rgmm(quantity)
        accepted_samples <- proposals[density(proposals) > runif(quantity, 0, gmm_bound(proposals))]
        samples <- c(samples, accepted_samples)
        
        pseudo_m <- quantity / length(accepted_samples)
        print(paste("Acceptance Rate: ", 1 / pseudo_m))
    }
    return(samples[1:num_samples])
}