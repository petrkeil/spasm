################################################################################
#
# Comparison of the pairwise and matrixwise ISA approaches.
#
# Petr Keil
#
################################################################################

library(spasm)

# the two matrices displayed in Fig. 2

# segregation
m1 <- matrix(c(1,0,0,0,
               0,1,0,0,
               0,0,1,0,
               1,1,1,1), byrow=TRUE, nrow=4, ncol=4)

# attraction
m2 <- matrix(c(1,0,0,0,
               1,0,0,0,
               1,0,0,0,
               1,1,1,1), byrow=TRUE, nrow=4, ncol=4)

#
m3 <- matrix(c(1,0,0,0,
               1,1,0,0,
               1,0,1,0,
               1,0,0,1), byrow=TRUE, nrow=4, ncol=4)


# ISA perspective
Whittaker(m1)
Whittaker(m2)
Whittaker(m3)


C_jacc(m1)
C_jacc(m2)
mean(C_jacc(m1))
mean(C_jacc(m2))

# Beta diversity perspective
Whittaker(t(m1))
Whittaker(t(m2))
C_jacc(t(m1))
C_jacc(t(m2))
mean(C_jacc(t(m1)))
mean(C_jacc(t(m2)))

# inverse of proportional fill
1/(sum(m1)/(nrow(m1)*ncol(m1)))


# mean numbers of species for SAR
colSums(m1)
mean(colSums(m1))
colSums(m2)
mean(colSums(m2))
