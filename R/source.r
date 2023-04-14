#' @export
mkData <- function(N1, N2, nitem, mu_2, var_2, uni_dif, nonuni_dif, anchor_n) {

  ipars <- genIRTpar(nitem)
  theta_1 <- rnorm(N1, 0, 1)

  # Make IRT data -----------------------------------------------------------
  # First group
  data_g1 <- simData.irt(ipar = ipars, theta = theta_1)
  data_g1 <- data.frame(data_g1)
  data_g1$group <- 0

  # Second group
  ipars$d <- ipars$d + uni_dif
  ipars$a <- ipars$a + nonuni_dif
  theta_2 <- rnorm(N2, mu_2, sqrt(var_2))

  data_g2 <- simData.irt(ipar = ipars, theta = theta_2)
  data_g2 <- data.frame(data_g2)
  data_g2$group <- 1

  data_irt <- rbind(data_g1, data_g2)

  # Make Stan data -----------------------------------------
  data_stan <- mk_staninpdata(data_irt, anchor_n = anchor_n)

  list(ipars = ipars,
       response = data_stan$lv.resp,
       group = data_stan$group,
       stan_data = data_stan[-1])
}


# select files
str_select <- function(out, regex, exclude = F) {
  if(exclude) {
    out[-grep(regex, out)] # str_detect
  } else {
    out[grep(regex, out)]
  }
}

# make index for DIF effects
get_difeffect_idx <- function(idx) {
  temp1 <- which(idx > 0)
  temp_length <- 1:length(temp1)

  difeffect_idx <- idx
  difeffect_idx[difeffect_idx > 0] <- temp1- temp_length
  difeffect_idx
}

# generate IRT data based on given theta
simData.irt <- function(ipar, theta) {

  a <- ipar$a
  d <- ipar$d
  g <- ipar$g

  nexaminee <- length(theta)
  nitem <- length(a)

  pr <- matrix(NA, nexaminee, nitem)
  for (j in 1:nitem){

    z = d[j]  + a[j]*theta

    pr[,j] <- 1 / (1 + exp(-z))

  }

  resp <- apply(pr, 2, function(x) {rbinom(nexaminee, 1, x)} )
  colnames(resp) <- paste0("u", 1:nitem)
  resp
}

# Convert matrix-like IRT data into Stan-like data.
mk_stantype <- function(itemdata) { # info = sim_info

  N = nrow(itemdata)
  nsec = ncol(itemdata)

  lv.resp <- itemdata

  total_N <- N
  # total_N <- N

  nworked <- rep(floor(nsec * 1), total_N)

  studentM <- do.call("c", lapply(seq(total_N),
                                  function(n) rep(n,each=nworked[n])))
  section <- do.call("c", lapply(seq(total_N),
                                 function(n) {
                                   sort(sample(1:nsec, nworked[n],
                                               replace = FALSE))}))
  ss <- cbind(studentM, section)
  grad <- sapply(1:dim(ss)[1], function(n) lv.resp[ss[n,1], ss[n,2]] )


  missing_pos <- is.na(grad)
  grad     <- grad[!missing_pos]
  studentM <- studentM[!missing_pos]
  section  <- section[!missing_pos]

  res <- list(
    lv.resp  = lv.resp,
    response = grad,       # score
    studIdx  = studentM,   # ppindex
    itemIdx  = section     # itemindex
  )

  return(res)
}

# Make Stan input data
mk_staninpdata <- function(dat, anchor_n, anchor_pos = 0) {

  if(anchor_pos[1] != 0) {
    if(anchor_n != length(anchor_pos)) {
      stop("anchor number not equal to anchor position!")
    }
  }

  nitem = ncol(dat)-1

  data_stan <- mk_stantype(dat[, -ncol(dat)])

  data_stan$group <- dat[, ncol(dat)]             # group
  data_stan$firstitem <- c(1, rep(0, nitem-1))    # firstitem
  data_stan$nitem <- nitem                         # nitem
  data_stan$nstud <- nrow(dat)                    # npp
  data_stan$nitemWorked <- length(data_stan$response)  # nppItem

  data_stan$unidif_n <- nitem - anchor_n          # unidif_n
  data_stan$nondif_n <- nitem - anchor_n          # nondif_n

  if(anchor_pos[1] == 0) {
    # unidif_idx
    data_stan$unidif_idx <- c(rep(0,anchor_n), rep(1, data_stan$unidif_n))
    data_stan$unidifeffect_idx <- get_difeffect_idx(data_stan$unidif_idx)

    # nondif_idx
    data_stan$nondif_idx <- c(rep(0,anchor_n), rep(1, data_stan$nondif_n))
    data_stan$nondifeffect_idx <- get_difeffect_idx(data_stan$nondif_idx)
  } else {

    # anchor_pos = c(1,2,3,6)

    data_stan$unidif_idx <- rep(1, nitem)
    data_stan$unidif_idx[anchor_pos] <- 0
    data_stan$nondif_idx <- rep(1, nitem)
    data_stan$nondif_idx[anchor_pos] <- 0

    # unidif_idx
    data_stan$unidifeffect_idx <- get_difeffect_idx(data_stan$unidif_idx)

    # nondif_idx
    data_stan$nondifeffect_idx <- get_difeffect_idx(data_stan$nondif_idx)
  }
  data_stan
}

# Cummulative product by row
rowprod <- function(x) {
  apply(x, 1, function(xx) Reduce(`*`, xx))
}

# Generate Theta
genTheta <- function(N, MU, SIG, BETA, sigE=1) {

  X <- mvrnorm(N, MU, SIG)
  X[,1] <- rbinom(N, 1, 0.5)

  ETA <- X %*% BETA + rnorm(N, 0, sigE)
  data <- data.frame(theta = ETA, X)
  data
}

# Generate IRT parameters
genIRTpar <- function(nitem=10) {

  a = round(rlnorm(nitem, 0.1, 0.5),3)
  a[1] <- 1
  a[a > 1.5 ] <- round(runif(length(a[a > 1.5 ]), 0.8, 1.2),3)
  a[a < 0.6 ] <- round(runif(length(a[a < 0.6 ]), 0.8, 1.2),3)
  # a
  d = round(rnorm(nitem, 0, 1),3)
  d[1] <- 0
  d[d > 2 ] <- round(runif(length(d[d > 2 ]), 1, 2),3)
  d[d < -2 ] <- round(runif(length(d[d < -2 ]), -2, -1),3)
  d[1] <- 0
  # d
  g = 0
  ipar <- data.frame(a, d, g)

  return(ipar)
}

# method for all dichotomous IRT model
simData.mg <- function(ipar, data, DIF_eff, D = 1){

  a <- ipar$a
  d <- ipar$d
  g <- ipar$g

  theta <- data[  grep("ID|theta|X", names(data))  ]
  X <- data[  grep("ID|X", names(data))  ]

  nitem <- length(a)

  theta0 <- theta[theta$X == 0, ]
  theta1 <- theta[theta$X == 1, ]

  # for x = 0
  nexaminee <- dim(theta0)[1]
  ID0 <- theta0$ID
  theta <- theta0$theta
  pr <- matrix(NA, nexaminee, nitem)
  for (j in 1:nitem){
    pr[,j] <- g[j] + (1 - g[j]) / (1 + exp(-D * (d[j] + a[j] * theta)))
  }
  resp <- apply(pr, 2, function(x) {rbinom(nexaminee, 1, x)} )
  resp0 <- cbind(ID0, resp)

  # for x = 1
  nexaminee <- dim(theta1)[1]
  ID1 <- theta1$ID
  theta <- theta1$theta
  pr <- matrix(NA, nexaminee, nitem)
  for (j in 1:nitem){
    pr[,j] <- g[j] + (1 - g[j]) / (1 + exp(-D * (d[j] + a[j] * theta)))
  }

  DIF_eff <- unlist(DIF_eff)
  for(ii in 1:length(DIF_eff)) {

    DIF_prod <- DIF_eff[ii]

    dif_item <- nitem - ii + 1

    a0 <- a[dif_item]
    d0 <- d[dif_item]
    g0 <- g[dif_item]

    pr[,dif_item] <- g0 + (1 - g0) / (1 + exp(-D * (DIF_prod + d0 + a0 * theta)))
  }
  resp <- apply(pr, 2, function(x) {rbinom(nexaminee, 1, x)} )
  resp1 <- cbind(ID1, resp)

  resp <- rbind(resp0,resp1)
  colnames(resp) <- c("ID", paste0("u", 1:nitem))
  resp <- resp[order(resp[,1]), ]

  data.frame(data, resp)
}

# method for all dichotomous IRT model
simData.dich <- function(ipar, data, DIF_eff, D = 1){

  a <- ipar$a
  d <- ipar$d
  g <- ipar$g

  theta <- data$theta
  X <- data[grep("X", names(data))]

  nexaminee <- length(theta)
  nitem <- length(a)

  pr <- matrix(NA, nexaminee, nitem)
  for (j in 1:nitem){
    pr[,j] <- g[j] + (1 - g[j]) / (1 + exp(-D * (d[j] + a[j] * theta)))
  }


  # uniform-DIF ---------------------------------------------------
  for(ii in 1:length(DIF_eff)) {

    DIF_prod <- as.matrix.data.frame(X) %*% DIF_eff[[ii]]

    dif_item <- nitem - ii + 1

    a0 <- a[dif_item]
    d0 <- d[dif_item]
    g0 <- g[dif_item]

    pr[,dif_item] <- g0 + (1 - g0) / (1 + exp(-D * (DIF_prod + d0 + a0 * theta)))
  }

  resp <- apply(pr, 2, function(x) {rbinom(nexaminee, 1, x)} )

  colnames(resp) <- paste0("u", 1:nitem)
  data.frame(data, resp)
}

