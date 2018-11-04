# 回転を定めるクオータニオンq1
library(onion)
theta <- runif(1) * 2*pi
ijk <- rnorm(3)
ijk <- ijk/sqrt(sum(ijk^2))

q1 <- cos(theta/2) + sin(theta/2) * (Hi*ijk[1]+Hj*ijk[2]+Hk*ijk[3])

# 別の回転 
theta <- runif(1) * 2*pi
ijk <- rnorm(3)
ijk <- ijk/sqrt(sum(ijk^2))

q2 <- cos(theta/2) + sin(theta/2) * (Hi*ijk[1]+Hj*ijk[2]+Hk*ijk[3])

# q1の回転をしてq2の回転をするのはq3
q3 <- q2 * q1
q3. <- Conj(q1) * Conj(q2)
Conj(q3)
q3.

# 回転q1と回転q3の「違い」は
q3/q1
q1/q3

# このq1/q3「違い」も回転で、それには、角度がある、それは
acos(Re(q3/q1))*2
# piを単位にすれば
acos(Re(q3/q1))*2/pi
