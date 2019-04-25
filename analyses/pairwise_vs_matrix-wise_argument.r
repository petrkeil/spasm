library(spasm)
?spasm

m1 <- matrix(c(1,0,0,0,
               0,1,0,0,
               0,0,1,0,
               1,1,1,1), byrow=TRUE, nrow=4, ncol=4)

m2 <- matrix(c(1,0,0,0,
               1,0,0,0,
               1,0,0,0,
               1,1,1,1), byrow=TRUE, nrow=4, ncol=4)

# ISA metrics
C_w(m1)
C_w(m2)
mean(C_jacc(m1))
mean(C_jacc(m2))

# Beta diversity metrics
C_w(t(m1))
C_w(t(m2))
mean(C_jacc(t(m1)))
mean(C_jacc(t(m2)))

# mean numbers of species for SAR
colSums(m1)
mean(colSums(m1))
colSums(m2)
mean(colSums(m2))
