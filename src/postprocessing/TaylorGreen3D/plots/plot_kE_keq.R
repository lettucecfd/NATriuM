
# require tikzDevice package to be installed
require (tikzDevice)
# load tikzDevice package  [if not available install.packages('...')]
library(tikzDevice)
# specify file name and size
tikz("dEdt_keq.tex", standAlone=TRUE, height=2.8, width=6.3)

# Reynolds number
Re=1600
nu=1./Re 
# read data
sim = read.table(file="../global_turbulence_table.txt")
# reference data
ref2 = read.table(file="Brachet_1983_1600.txt")
# to create derivatives.txt, call 'python calculate_dEdt_inc.py'
sim2 = read.table(file="derivatives.txt")
# turbulent dissipation
sim$V3 = nu / (2*pi)**3 *  (sim$V54 + sim$V63 + sim$V73 + sim$V84 + sim$V96 + sim$V109 + sim$V123 + sim$V138+ sim$V154)
# enstrophy
sim$V4 = (sim$V138 - 2* sim$V136 + sim$V109 + sim$V73 - 2*sim$V119 + sim$V123 + sim$V84 - 2* sim$V82 + sim$V63) / (2*pi)**3

# increment in the finite differences (has to agree with 'n' in calculate_dEdt_inc.py)
N = 5
# overwrite sim$V5 by finite differences
A = -sim2$V3
B = array(NA,N+1)
sim$V5=append(A,B)


# margin (bottom, left, top, right)
par( mai = c(0.7, 0.8, 0.2, 2.5), xpd=TRUE )


# plot
plot(sim$V2, sim$V3, xlab="Time $t$", ylab="Dissipation", type="l", lwd=4, lty=1, ylim=c(0.0, 1.2*max(sim$V5,na.rm=TRUE)))
points(ref2, type ="l", lwd=2, col="red", pch=10, cex=1.7)
points(sim$V2, sim$V5, type ="l", lwd=2, col="blue", lty=2)
points(sim$V2[seq(1,length(sim$V2),length(sim$V2)/14)], nu * sim$V4[seq(1,length(sim$V2),length(sim$V2)/14)], type ="p", lwd=1, col="cyan2", lty=1, pch=15, cex=0.8)
points(sim$V2, nu * sim$V4, type ="l", lwd=1, col="cyan2", lty=1, pch=15)
points(sim$V2, sim$V5-sim$V3, type ="l", lwd=2, col="orange", lty=3)

legend("topright", inset=c(-0.75,0), legend = c("Kin. Energy Dissipation $-\\frac{dk}{dt}$", "Turbulent Dissipation $\\epsilon$", "Scaled Enstrophy $\\nu \\mathcal{E}$", "Numerical Dissipation","Reference"), pch=c(-1,-1,15,-1,-1), lty=c(2,1,1,3,1), lwd=c(2,4,1,2,2), col=c("blue","black","cyan2","orange","red"), pt.cex=c(1,1,0.8,1,1), bty="n")


# write data to files for quicker processing
R = sim$V4/sim$V5
eff_Re = mean(R[sim$V2<20],na.rm=TRUE)
write(eff_Re, "eff_Re.txt")
min_Re = min(sim$V4/sim$V5,na.rm=TRUE)
write(min_Re, "min_Re.txt")

data_write = data.frame(sim$V2, sim$V4, sim$V5)
colnames(data_write) = c("t", "enstrophy", "dissipation")
write.table(data_write, file="data.txt")
# finalize
dev.off()


