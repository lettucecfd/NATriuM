require (tikzDevice)
library(tikzDevice)
tikz("eff_Re.tex", standAlone=TRUE, height=2.8, width=2.8)

#source("log_axis.R")

Re = c(100, 200, 400, 800, 1600, 3000, 5000, 15000, 100000, 1e6)
bgk = c(100.4955, 198.5857, 389.9639, 740.4227, NA, NA, NA, NA, NA, NA)
present = c(100.4155, 197.3707, 375.8985, 663.9783, 1044.21, 1433.091, 1774.389, 2514.173, 3679.024, 4205.941)

#ref2 = read.table(file="Brachet_1983_800.txt") #just to make sure Gassner and Beck extracted the right data from the plot

# margin (bottom, left, top, right)
par( mai = c(0.7, 0.8, 0.2, 0.2) )
source("colors.R")
# plot

plot (1, type="n", xaxt="n", yaxt="n", xlab = "Reynolds Number $Re$", ylab = "Effective Reynolds Number $\\langle Re_\\mathrm{eff} \\rangle_t$", xlim=c(80,20000), ylim=c(80,20000), log="xy")
#minor.ticks.axis(1,3,mn=1,mx=5)
axis(1,at=c(100,1000,10000), labels=c("100","1000","10000"))
axis(1,at=c(seq(100,1000,100),seq(1000,10000,1000),seq(10000,100000,10000)),labels=FALSE,tcl=-0.25)

axis(2,at=c(100,1000,10000), labels=c("100","1000","10000"))
axis(2,at=c(seq(100,1000,100),seq(1000,10000,1000),seq(10000,100000,10000)),labels=FALSE,tcl=-0.25)

abline(0,1, col="black", lwd=1)
points (Re,bgk, lwd=3, lty=1, pch=3, type="o", col=BGKC1, cex=2)
points(Re,present, lwd=3, lty=1, pch=1, type="o", col=REGC1)
legend("topleft", legend=c("BGK", "Stabilized"), pch = c(3,1), lty=c(1,1), lwd=c(3,3), col=c(BGKC1,REGC1), pt.cex=c(2,1), bty ="n")
#plot(sim$V2, sim$V3, xlab="Time $t$", ylab="Dissipation", type="l", lwd=4, lty=1, ylim=c(0.0, 1.2*max(ref2$V2)))
#points(ref2, type ="p", lwd=2, col="red", pch=10, cex=1.7)
#points(sim$V2, sim$V5, type ="l", lwd=2, col="blue", lty=2)
#points(sim$V2[seq(1,length(sim$V2),length(sim$V2)/14)], nu * sim$V4[seq(1,length(sim$V2),length(sim$V2)/14)], type ="p", lwd=1, col="cyan2", lty=1, pch=15, cex=0.8)
#points(sim$V2, nu * sim$V4, type ="l", lwd=1, col="cyan2", lty=1, pch=15)
#points(sim$V2, sim$V5-sim$V3, type ="l", lwd=2, col="orange", lty=3)
#legend(x=min(sim$V2), y=1.2*max(sim$V3), legend = c("Kin. Energy Dissipation $-\\frac{dk}{dt}$", "Turbulent Dissipation $\\epsilon$", "Scaled Enstrophy $\\nu \\mathcal{E}$", "Numerical Dissipation","Re
#ference"), pch=c(-1,-1,15,-1,10), lty=c(2,1,1,3,-1), lwd=c(2,4,1,2,2), col=c("blue","black","cyan2","orange","red"), pt.cex=c(1,1,0.8,1,1.7), bty="n")

dev.off()


