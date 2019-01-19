

## kernel density function
kd_fun <- function(i, dat, md){
  # treatment
  den.grp <- density(dat[,i])
  datfram.grp <- data.frame(x = den.grp$x, y = den.grp$y)
  ix = which.min(abs(den.grp$x - md[i]))
  datfram.grp[ix,"y"]
}
