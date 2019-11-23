ID <-  c(1,6,3,0,0,2,2)
density <- function(x) {
    t1 <- ID[1] * exp(-sin((ID[1] * x ** 2)/(15 - ID[1]) - (x - 3 - ID[2] * pi) ** 2 / (2 * (5 + ID[3]) ** 2)))
    t2 <- 2 * (1 + ID[7]) * exp(- x**2 /32)
    t3 <- (10 - ID[7]) * exp(-cos(ID[4] * x ** 2 / (15 + ID[4])) - (x + 3 + ID[5] * pi)**2/(2*(5+ID[6])**2))
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