# --------------------------------------------------------------------------------
# Useful functions used in this project
# --------------------------------------------------------------------------------

### Life table based on mx by single ages and by sex

AKm02a0 <- function(m0, sex = "m"){
  sex <- rep(sex, length(m0))
  ifelse(sex == "m",
         ifelse(m0 < .0230, {0.14929 - 1.99545 * m0},
                ifelse(m0 < 0.08307, {0.02832 + 3.26201 * m0},.29915)),
         # f
         ifelse(m0 < 0.01724, {0.14903 - 2.05527 * m0},
                ifelse(m0 < 0.06891, {0.04667 + 3.88089 * m0}, 0.31411))
  )
}

life.table <- compiler::cmpfun(function(mx,sex = "f"){
  i.openage <- length(mx)
  OPENAGE   <- i.openage - 1
  RADIX     <- 1
  ax        <- mx * 0 + .5
  ax[1]     <- AKm02a0(m0 = mx[1], sex = sex)
  qx        <- mx / (1 + (1 - ax) * mx)
  qx[i.openage]       <- ifelse(is.na(qx[i.openage]), NA, 1)
  ax[i.openage]       <- 1 / mx[i.openage]
  px 				    <- 1 - qx
  px[is.nan(px)]      <- 0
  lx 			        <- c(RADIX, RADIX * cumprod(px[1:OPENAGE]))
  dx 				    <- lx * qx
  Lx 				    <- lx - (1 - ax) * dx
  Lx[i.openage ]	    <- lx[i.openage ] * ax[i.openage ]
  Tx 				    <- c(rev(cumsum(rev(Lx[1:OPENAGE]))),0) + Lx[i.openage]
  ex 				    <- Tx / lx
  #ex[1] uncomment this if you want only life expectancy at birth
  return(data.frame(qx=qx, px = px, ax = ax, lx = lx , dx = dx, Lx= Lx,
                    Tx = Tx, ex = ex))
})


# edag = ineq_edag(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax)
# life disparity

ineq_edag <- function(age, dx, lx, ex, ax){
  age_length_equal <- all.equal(length(age),length(dx),
                                length(lx),length(ex),
                                length(ax))

  stopifnot(age_length_equal)

  # length of the age interval
  n <- c(diff(age),1)
  explusone <- c(ex[-1],ex[length(age)])
  # the average remaining life expectancy in each age interval
  # (as opposed to the beginning of the interval)
  # ends up being roughly half of the ex between ages
  ex_average <- ex + ax / n * (explusone - ex)

  rev(cumsum(rev(dx * ex_average))) / lx
}


# for more than one time being saved:

ineq_edag_i <- function(age, dx, lx, ex, ax,i){
  age_length_equal <- all.equal(length(age),length(dx),
                                length(lx),length(ex),
                                length(ax))

  stopifnot(age_length_equal)

  # length of the age interval
  n <- c(diff(age),1)
  explusone <- c(ex[-1],ex[length(age)])
  # the average remaining life expectancy in each age interval
  # (as opposed to the beginning of the interval)
  # ends up being roughly half of the ex between ages
  ex_average <- ex + ax / n * (explusone - ex)

  (rev(cumsum(rev(dx * ex_average*(-log(lx)^(i-1)))))/factorial(i))/lx
}



# Calculate a lifetable entropy H


ineq_H <- function(age, dx, lx, ex, ax){
  ineq_edag(age, dx, lx, ex, ax) / (ex + age)
}

# adapting for revivorship

ineq_H_i <- function(age, dx, lx, ex, ax,i){
  ineq_edag_i(age, dx, lx, ex, ax,i) / (ex + age)
}







