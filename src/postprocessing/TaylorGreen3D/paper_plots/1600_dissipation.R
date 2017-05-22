require (tikzDevice)
library(tikzDevice)
tikz("1600_dissipation.tex", standAlone=TRUE, height=3.2, width=4)

#set palette
source("colors.R")
CLR = tol15rainbow
PTS = pch4lines

#read data
Re100_bgk = read.table(file="../results/Re100-ref4-p5-coll0-sl1-CFL7.07107/plot/data.txt")
Re100_reg = read.table(file="../results/Re100-ref4-p5-coll0-sl1-CFL7.07107-reg1/plot/data.txt")
Re200_bgk = read.table(file="../results/Re200-ref4-p5-coll0-sl1-CFL7.07107/plot/data.txt")
Re200_reg = read.table(file="../results/Re200-ref4-p5-coll0-sl1-CFL7.07107-reg1/plot/data.txt")
Re400_bgk = read.table(file="../results/Re400-ref4-p5-coll0-sl1-CFL7.07107/plot/data.txt")
Re400_reg = read.table(file="../results/Re400-ref4-p5-coll0-sl1-CFL7.07107-reg1/plot/data.txt")
Re800_bgk = read.table(file="../results/Re800-ref4-p5-coll0-sl1-CFL7.07107/plot/data.txt")
Re800_reg = read.table(file="../results/Re800-ref4-p5-coll0-sl1-CFL7.07107-reg1/plot/data.txt")
Re1600_bgk = read.table(file="../results/Re1600-ref4-p5-coll0-sl1-CFL7.07107/plot/data.txt")
Re1600_reg = read.table(file="../results/Re1600-ref4-p5-coll0-sl1-CFL7.07107-reg1/plot/data.txt")

Re100_ref = read.table(file="Brachet_1983_100.txt")
Re200_ref = read.table(file="Brachet_1983_200.txt")
Re400_ref = read.table(file="Brachet_1983_400.txt")
Re800_ref = read.table(file="Brachet_1983_800.txt")
Re1600_ref = read.table(file="Brachet_1983_1600.txt")

#preparation
seq_14=seq(1,length(Re100_bgk$t),length(Re100_bgk$t)/14)
seq_12=seq(1,length(Re100_bgk$t),length(Re100_bgk$t)/12)

#plot (bootom left top right)
par( mai = c(0.7, 0.8, 0.2, 0.2) )

plot(1,type="n", xlim=c(0,20), ylim=c(0,0.014), xlab="Time $t$", ylab="Dissipation")
#Re100
#points (Re100_bgk$t, Re100_bgk$dissipation, type="l", lw=1, col="dodgerblue4")
#points (Re100_bgk$t[seq_14], Re100_bgk$dissipation[seq_14], type="p", lw=3, pch=3, col="dodgerblue4")
#points (Re100_reg$t, Re100_reg$dissipation, type="l", lw=1, col="darkgoldenrod2")
#points (Re100_reg$t[seq_12], Re100_reg$dissipation[seq_12], type="p", lw=3, pch=1, col="darkgoldenrod2")
#points (Re100_ref$V1, Re100_ref$V2, type="l", lw=3, col="black")
#Re200
#points (Re200_bgk$t, Re200_bgk$dissipation, type="l", lw=1, col="dodgerblue4")
#points (Re200_bgk$t[seq_14], Re200_bgk$dissipation[seq_14], type="p", lw=3, pch=3, col="dodgerblue4")
#points (Re200_reg$t, Re200_reg$dissipation, type="l", lw=1, col="darkgoldenrod2")
#points (Re200_reg$t[seq_12], Re200_reg$dissipation[seq_12], type="p", lw=3, pch=1, col="darkgoldenrod2")
#points (Re200_ref$V1, Re200_ref$V2, type="l", lw=3, col="black")
#Re400
#points (Re400_bgk$t, Re400_bgk$dissipation, type="l", lw=1, col="dodgerblue4")
#points (Re400_bgk$t[seq_14], Re400_bgk$dissipation[seq_14], type="p", lw=3, pch=3, col="dodgerblue4")
#points (Re400_reg$t, Re400_reg$dissipation, type="l", lw=1, col="darkgoldenrod2")
#points (Re400_reg$t[seq_12], Re400_reg$dissipation[seq_12], type="p", lw=3, pch=1, col="darkgoldenrod2")
#points (Re400_ref$V1, Re400_ref$V2, type="l", lw=3, col="black")
#Re800
points (Re1600_bgk$t, Re1600_bgk$dissipation, type="l", lw=2, col=BGKC1)
points (Re1600_bgk$t[seq_14], Re1600_bgk$dissipation[seq_14], type="p", lw=2, pch=BGKP1, col=BGKC1)
points (Re1600_reg$t, Re1600_reg$dissipation, type="l", lw=2, col=REGC1)
points (Re1600_reg$t[seq_12], Re1600_reg$dissipation[seq_12], type="p", lw=2, pch=REGP1, col=REGC1)
points (Re1600_bgk$t, Re1600_bgk$enstrophy/1600, type="l", lw=2, col=BGKC2)
points (Re1600_bgk$t[seq_14], Re1600_bgk$enstrophy[seq_14]/1600, type="p", lw=2, pch=BGKP2, col=BGKC2)
points (Re1600_reg$t, Re1600_reg$enstrophy/1600, type="l", lw=2, col=REGC2)
points (Re1600_reg$t[seq_12], Re1600_reg$enstrophy[seq_12]/1600, type="p", lw=2, pch=REGP2, col=REGC2)
points (Re1600_ref$V1, Re1600_ref$V2, type="l", lw=4, lty = 2, col="black")
#Re1600
#points (Re1600_bgk$t, Re1600_bgk$dissipation, type="l", lw=1, col="dodgerblue4")
#points (Re1600_bgk$t[seq_14], Re1600_bgk$dissipation[seq_14], type="p", lw=3, pch=3, col="dodgerblue4")
#points (Re1600_reg$t, Re1600_reg$dissipation, type="l", lw=1, col="darkgoldenrod2")
#points (Re1600_reg$t[seq_12], Re1600_reg$dissipation[seq_12], type="p", lw=3, pch=1, col="darkgoldenrod2")
#points (Re1600_ref$V1, Re1600_ref$V2, type="l", lw=3, col="black")

#legend("topright", bty="n", legend=c("Stabilized $-dk/dt$","Stabilized $\\nu \\mathcal{E}$","BGK $-dk/dt$","BGK $\\nu \\mathcal{E}$","Reference"), col=c( "darkgoldenrod2","darkgoldenrod2","dodgerblue4","dodgerblue4", "black"), lw=c(3,3,3,3,3), lty=c(-1,-1,-1,-1,1), pch=c(1,10,3,17,-1))
title("$Re=1600$")
#text(-0.5,0.004,"$Re=800$")
dev.off()
