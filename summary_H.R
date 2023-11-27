setwd('/Users/jose/Documents/Proyectos/Entropy/spike14_manual/Statistics')

install.packages("plot.matrix")
library("plot.matrix")

rm(list = ls())

#-------------------------------------------------------------------------------

dset1 = read.table('statistics_1.tsv', header=TRUE, sep='\t')
summary(dset1)

plot(dset1$Ipct ~ dset1$fragment, lty=1, lwd=2, cex=0.4, col='black', type = 'o', xlab='alignment position', ylab='%', ylim=c(0,100), main='Analysis of SARS-CoV-2')
points(dset1$Dpct ~ dset1$fragment, lty=1, lwd=2, cex=0.4, col='red', type = 'o')
points((dset1$coverage)*100 ~ dset1$fragment, lty=1, lwd=2, cex=0.4, col='green', type = 'o')
legend(1, 80, legend = c('Ipct', 'Dpct', 'coverage'), col = c("black","red","green"), lty = 1, cex = 0.7, lwd = 2, bg = 'white')

dset2 = read.table('statistics_2.tsv', header=TRUE, sep='\t')
summary(dset2)

unique(dset2$frgmt_1)
length(unique(dset2$frgmt_1))
unique(dset2$frgmt_2)
unique(dset1$fragment)
m <- matrix(0, nrow = length(unique(dset2$frgmt_1)), ncol = length(unique(dset2$frgmt_1)), dimnames = list((c(unique(dset2$frgmt_1))),c(unique(dset2$frgmt_1))))
#m <- matrix(0, nrow = 60, ncol = 60) # Matriz de prueba
#runif(1) # Genera numeros aleatorios entre 0 y 1
colnames(m)
rownames(m)
nrow(m)
ncol(m)

vec <- unique(dset2$frgmt_1)
x = 1
y = 1
r = 0
for (i in 1:length(unique(dset2$frgmt_1))){
  for (j in i:length(unique(dset2$frgmt_1))){
    #print (i)
    #print (j)
    ##print (vec[i])
    ##print (vec[j])
    r = r + 1
    #print (dset2$Ipct[r])
    #m[i,j] <- dset2$Ipct[r]
    #print (dset2$Dpct[r])
    m[i,j] <- dset2$Dpct[r]
    #print (dset2$unique[r])
    #m[i,j] <- dset2$unique[r]
    #print (dset2$coverage[r])
    #m[i,j] <- dset2$coverage[r]
    ##m[i,j] <- runif(1)
  }
}
m
summary(m)
m[,1] # columnas
m[1,] # renglones
heatmap(m)
#pdf("statistics_2_D_250_200.pdf")
par(mar=c(5.1, 4.5, 4.5, 4.5))
plot(m, col = rainbow(10, rev = TRUE, start = c(0,0.2), end = c(0.8,1)), cex.axis = 0.7, main = "Phylogenetic dissonance (D)", xlab = "column number along the alignment", ylab = "column number along the alignment")
plot(m, col = rainbow(10, rev = TRUE, start = c(0,0.2), end = c(0.8,1)), main = "Coverage", xlab = "column number along the alignment", ylab = "column number along the alignment")
plot(m, col = rainbow(10, rev = TRUE, start = c(0,0.2), end = c(0.8,1)), main = "Information", xlab = "column number along the alignment", ylab = "column number along the alignment")
plot(m, col = topo.colors)
# este es el bueno
plot(m, col = rainbow(15, rev = TRUE, start = c(0,1), end = c(0.5,1)), las = 2, cex = 0.3, main = "Phylogenetic dissonance (D)", xlab = "column number along the alignment", ylab = "column number along the alignment")
plot(m, col = rainbow(8, rev = TRUE, start = c(0,0.2), end = c(0.8,1)), las = 2, cex = 0.3)
#dev.off()
?rainbow()
?par

m[57,1] <- 0
m[56,1] <- 11
m[55,1] <- 24
m[54,1] <- 37
m[53,1] <- 49
m[52,1] <- 61
m[51,1] <- 74
m[50,1] <- 87
m[49,1] <- 100

((1596 * (60+25))/60)/24
(1596 * (60+25)) # numero de segundos totales
(1596 * (60+25))/60 # numero de minutos totales 
((1596 * (60+25))/60)/60 # numero de horas totales
(((1596 * (60+25))/60)/60)/24 # numero de dias totales

# Cuantas celdas
(58^2/2) - 58/2 
1682 - 1653
(58^2/2) - 58
57+56+55

v <- c()
for (i in 1:57){
  print (i)
  v[i] <- i
}
v
sum(v)
1653*2 + 58 == 58^2

#-------------------------------------------------------------------------------

dsetss = read.table('statistics_ss_manual.tsv', header=TRUE, sep='\t')
summary(dsetss)
dsetss$lnM1 <- dsetss$lnL1 + dsetss$lnL2
# lnB10 = 2*(lnM1 - lnM0)
#dsetss$lnB10 <- 2*(dsetss$lnM1 - dsetss$lnM0) 
dsetss$lnB10 <- 2*(dsetss$lnM0 - dsetss$lnM1) # indicacion de Lizbeth
# K = lnM0 - lnM1
# if K >  1, then model M0 wins (CONCATENATED)
# if K < -1, then model M1 wins (SEPARATED)
dsetss$K <- dsetss$lnM0 - dsetss$lnM1
# Criteria (2) from Kass and Raftery (1995)
# exp(lnM1 - lnM0)
dsetss$KR <- exp(dsetss$lnM1 - dsetss$lnM0)

unique(dsetss$frgmt_1s)
length(unique(dsetss$frgmt_1s))
unique(dsetss$frgmt_2s)
length(unique(dsetss$frgmt_2s))
mss <- matrix(0, nrow = length(unique(dsetss$frgmt_1s)), ncol = length(unique(dsetss$frgmt_1s)), dimnames = list((c(unique(dsetss$frgmt_1s))),c(unique(dsetss$frgmt_1s))))
colnames(mss)
rownames(mss)
nrow(mss)
ncol(mss)

vec <- unique(dsetss$frgmt_1s)
x = 1
y = 1
r = 0
for (i in 1:length(unique(dsetss$frgmt_1s))){
  for (j in i:length(unique(dsetss$frgmt_1s))){
    #print (i)
    #print (j)
    ##print (vec[i])
    ##print (vec[j])
    r = r + 1
    #print (dsetss$lnB10[r])
    #mss[i,j] <- dsetss$lnB10[r]
    #if (1 < 0){
    if (is.na(dsetss$lnB10[r])){
      mss[i,j] <- dsetss$lnB10[r]
    } else {
      if (dsetss$lnB10[r] > 10){
        mss[i,j] <- 12
      } else if (dsetss$lnB10[r] > 6 & dsetss$lnB10[r] <= 10){
        mss[i,j] <- 8
      } else if (dsetss$lnB10[r] > 2 & dsetss$lnB10[r] <= 6){
        mss[i,j] <- 4
      } else if (dsetss$lnB10[r] > 0 & dsetss$lnB10[r] <= 2){
        mss[i,j] <- 1
      } else if (dsetss$lnB10[r] == 0){
        mss[i,j] <- 0
      } else if (dsetss$lnB10[r] >= -2 & dsetss$lnB10[r] < 0){
        mss[i,j] <- -1
      } else if (dsetss$lnB10[r] >= -6 & dsetss$lnB10[r] < -2){
        mss[i,j] <- -4
      } else if (dsetss$lnB10[r] >= -10 & dsetss$lnB10[r] < -6){
        mss[i,j] <- -8
      } else if (dsetss$lnB10[r] >= -100 & dsetss$lnB10[r] < -10){
        mss[i,j] <- -50
      } else if (dsetss$lnB10[r] >= -1000 & dsetss$lnB10[r] < -100){
        mss[i,j] <- -500
      } else if (dsetss$lnB10[r] < -1000){
        mss[i,j] <- -5000
      }
    }
    #}
    #mss[i,j] <- dsetss$K[r]
    #mss[i,j] <- dsetss$KR[r]
  }
}
#for (i in 1:length(unique(dsetss$frgmt_1s))){
#  for (j in i:length(unique(dsetss$frgmt_1s))){
#    #print (dsetss$lnB10[r])
#    mss[j,i] <- 'NA'
#  }
#}
mss
summary(mss)
mss[,6] # columnas
mss[1,] # renglones
heatmap(mss)
#pdf("statistics_ss_250_200.pdf")
par(mar=c(5.1, 4.1, 4.1, 5.1))
plot(mss, cex.axis = 0.7, col = rainbow(25, rev = TRUE, start = c(0,0.2), end = c(0.8,1)), main = "Bayes Factors", xlab = "column number along the alignment", ylab = "column number along the alignment")
plot(mss, cex.axis = 0.7, main = "Bayes Factors", xlab = "column number along the alignment", ylab = "column number along the alignment")
plot(mss, col = topo.colors)
# Este es el bueno
plot(mss, col = rainbow(30, rev = TRUE, start = c(0,0.2), end = c(0.8,1)), las = 2, cex = 0.3, main = "Bayes Factors", xlab = "column number along the alignment", ylab = "column number along the alignment")
#dev.off()
?par
mss[57,1] <- 0
mss[56,1] <- 0
mss[55,1] <- 0
mss[54,1] <- 0
mss[53,1] <- 0
mss[52,1] <- 0
# -----
mss[50,1] <- -2000
mss[49,1] <- -1000
mss[48,1] <- -100
mss[47,1] <- -10
mss[46,1] <- 0
mss[45,1] <- 10
mss[44,1] <- 100
# -----
mss[43,1] <- 1000
mss[42,1] <- 5000
mss[41,1] <- 10000
