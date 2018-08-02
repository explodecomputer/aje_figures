#######################################################
# Simulation to illustrate the power                  #
# of Cochran's Q statistic, Rucker's Q statitistic    #
# and MR-Egger to detect global pleiotropy            #
#######################################################

rm(list=ls(all=TRUE))

set.seed(1)

# Simulation with no heterogeneity

QIVW1 = QEG   = NULL
Qp1   = Qp2   = Qp3 =  NULL

Power  = matrix(nrow=20,ncol=3)

N        = 5000
UB       = 0.095 #c(0.095,0.156,0.24,0.4,1)
b        = 0.05
L        = 25
U        = seq(0,0.2,length=20)

for(k in 1:20){
for(i in 1:N){

alphaJ1  = runif(L,-U[k],U[k])
L        = 25 ; DF = L-1
seBetaXG = runif(L,0.06,UB)
seBetaYG = runif(L,0.015,0.11)
BXG      = runif(L,0.34,1.1) 
BetaXG   = BXG + rnorm(L,0,seBetaXG)
Beta     = b
BetaYG   = Beta*BXG + alphaJ1 + rnorm(L,0,seBetaYG)
BIV      = BetaYG/BetaXG
W1       = 1/(seBetaYG^2/BetaXG^2)
BIVw1    = BIV*sqrt(W1)
sW1      = sqrt(W1)


# 1st order weights IVW

IVWfitR1  = summary(lm(BIVw1 ~ -1+sW1))
phi_IVW1  = IVWfitR1$sigma^2
QIVW1[i]  = DF*phi_IVW1
Qp1[i]    = 1-pchisq(QIVW1[i],DF)


# 1st order weights MR-Egger

EGfitR1   = summary(lm(BetaYG ~ BetaXG,weights=1/seBetaYG^2))
DF        = L-2
phi_IVW2  = EGfitR1$sigma^2
QEG[i]    = DF*phi_IVW2
Qp2[i]    = 1-pchisq(QEG[i],DF)
Qp3[i]    = EGfitR1$coef[1,4]
}

P1 = length(Qp1[Qp1<=0.05])/N
P2 = length(Qp2[Qp2<=0.05])/N
P3 = length(Qp3[Qp3<=0.05])/N

Power[k,] = c(P1,P2,P3)
print(k)

}

## Power

par(mfrow=c(1,2))


plot(U,Power[,1],ylim=c(0,1),
xlab="Maximum size of direct/pleiotropic effect",ylab="Power to detect global pleiotropy",
type="l",lty=1,lwd=3,cex.lab=1.5)
lines(U,Power[,2],lty=2,lwd=3)
lines(U,Power[,3],lty=3,lwd=3)
legend("right",c("Cochran's Q", "Rucker's Q", "MR-Egger intercept"),bty="n",
lty=c(1,2,3),lwd=3,cex=1.5)

lines(c(-1,0.2),rep(0.05,2),lty=4,lwd=1)

Figure2Adata = data.frame(U,Power[,1],Power[,2],Power[,3])

colnames(Figure2Adata) = c("MaxPleioSize","PowerCochrQ","PowerRuckQ","PowerMREgg")
Figure2Adata

write.csv(Figure2Adata, file = "Figure2Adata.csv")

x11() 

##################################################
# 2nd plot to show the individual contribution   #
# to Cochran's Q and Rucker's Q for a single     #
# simulation                                     #
##################################################


set.seed(2)
L        = 25 ; DF = L-1

alphaJ1   = runif(L,-U[5],U[5])

seBetaXG = runif(L,0.06,UB)
seBetaYG = runif(L,0.015,0.11)
BXG      = runif(L,0.34,1.1) 
BetaXG   = BXG + rnorm(L,0,seBetaXG)
Beta     = b
BetaYG   = Beta*BXG + 0.1 + alphaJ1 + rnorm(L,0,seBetaYG)
BIV      = BetaYG/BetaXG
W1       = 1/(seBetaYG^2/BetaXG^2)
BIVw1    = BIV*sqrt(W1)
sW1      = sqrt(W1)

# 1st order weights IVW

IVWfitR1  = summary(lm(BIVw1 ~ -1+sW1))
phi_IVW1  = IVWfitR1$sigma^2
QIVW1[i]  = DF*phi_IVW1
Qp1[i]    = 1-pchisq(QIVW1[i],DF)

# 1st order weights MR-Egger

EGfitR1   = summary(lm(BetaYG ~ BetaXG,weights=1/seBetaYG^2))
DF        = L-2
phi_IVW2  = EGfitR1$sigma^2
QEG[i]    = DF*phi_IVW2
Qp2[i]    = 1-pchisq(QEG[i],DF)
Qp3[i]    = EGfitR1$coef[1,4]

Q1ind     = W1*(BIV - IVWfitR1$coef[1])^2
Q2ind     = W1*(BIV - EGfitR1$coef[1,1]/BetaXG - EGfitR1$coef[2,1])^2

plot(Q1ind,pch=19,cex = 1.5,xlab="SNP",ylab="Q statistic contribution",
cex.lab=1.5,xlim=c(0,29))
points(Q2ind,pch=2,cex=1.5)
legend("topleft",bty="n",c("Cochran's Q (IVW)","Rucker's Q (MR-Egger)"),pch=c(19,2),cex=1.5)

lines(c(0,25),rep(qchisq(0.95,1),2),lty=2)
lines(c(0,25),rep(qchisq(1-(0.05/25),1),2),lty=3)
legend(24,4.8,"p=0.05",bty="n")
legend(24,10.3,"p=0.05/25",bty="n")


Figure2Bdata = data.frame(Q1ind,Q2ind)

colnames(Figure2Bdata) = c("CochranQ_contrib","RuckerQ_contrib")
Figure2Bdata

write.csv(Figure2Bdata, file = "Figure2Bdata.csv")

