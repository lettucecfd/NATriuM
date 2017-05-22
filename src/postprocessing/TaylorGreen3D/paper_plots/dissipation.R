require (tikzDevice)
library(tikzDevice)
tikz("dissipation.tex", standAlone=TRUE, height=3.2, width=4)


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
source("colors.R")

#plot (bootom left top right)
par( mai = c(0.7, 0.8, 0.2, 0.2) )

plot(1,type="n", xlim=c(-3,20), ylim=c(0,0.014), xlab="Time $t$", ylab="Dissipation $-dk/dt$")
#Re100
points (Re100_bgk$t, Re100_bgk$dissipation, type="l", lw=2, col=BGKC1)
points (Re100_bgk$t[seq_14], Re100_bgk$dissipation[seq_14], type="p", lw=2, pch=BGKP1, col=BGKC1)
points (Re100_reg$t, Re100_reg$dissipation, type="l", lw=2, col=REGC1)
points (Re100_reg$t[seq_12], Re100_reg$dissipation[seq_12], type="p", lw=2, pch=REGP1, col=REGC1)
points (Re100_ref$V1, Re100_ref$V2, type="l", lw=4, lty=2, col="black")
#Re200
points (Re200_bgk$t, Re200_bgk$dissipation, type="l", lw=2, col=BGKC1)
points (Re200_bgk$t[seq_14], Re200_bgk$dissipation[seq_14], type="p", lw=2, pch=BGKP1, col=BGKC1)
points (Re200_reg$t, Re200_reg$dissipation, type="l", lw=2, col=REGC1)
points (Re200_reg$t[seq_12], Re200_reg$dissipation[seq_12], type="p", lw=2, pch=REGP1, col=REGC1)
points (Re200_ref$V1, Re200_ref$V2, type="l", lw=4, lty=2, col="black")
#Re400
points (Re400_bgk$t, Re400_bgk$dissipation, type="l", lw=2, col=BGKC1)
points (Re400_bgk$t[seq_14], Re400_bgk$dissipation[seq_14], type="p", lw=2, pch=BGKP1, col=BGKC1)
points (Re400_reg$t, Re400_reg$dissipation, type="l", lw=2, col=REGC1)
points (Re400_reg$t[seq_12], Re400_reg$dissipation[seq_12], type="p", lw=2, pch=REGP1, col=REGC1)
points (Re400_ref$V1, Re400_ref$V2, type="l", lw=4, lty=2, col="black")
#Re800
#points (Re800_bgk$t, Re800_bgk$dissipation, type="l", lw=2, col=BGKC1)
#points (Re800_bgk$t[seq_14], Re800_bgk$dissipation[seq_14], type="p", lw=2, pch=BGKP1, col=BGKC1)
#points (Re800_reg$t, Re800_reg$dissipation, type="l", lw=2, col=REGC1)
#points (Re800_reg$t[seq_12], Re800_reg$dissipation[seq_12], type="p", lw=2, pch=REGP1, col=REGC1)
#points (Re800_ref$V1, Re800_ref$V2, type="l", lw=2, col="black")
#Re1600
#points (Re1600_bgk$t, Re1600_bgk$dissipation, type="l", lw=2, col=BGKC1)
#points (Re1600_bgk$t[seq_14], Re1600_bgk$dissipation[seq_14], type="p", lw=2, pch=BGKP1, col=BGKC1)
#points (Re1600_reg$t, Re1600_reg$dissipation, type="l", lw=2, col=REGC1)
#points (Re1600_reg$t[seq_12], Re1600_reg$dissipation[seq_12], type="p", lw=2, pch=REGP1, col=REGC1)
#points (Re1600_ref$V1, Re1600_ref$V2, type="l", lw=2, col="black")

legend("topright", bty="n", legend=c("BGK","Stabilized","Reference"), col=c(BGKC1, REGC1, "black"), lw=c(2,2,4), lty=c(1,1,2), pch=c(BGKP1,REGP1,-1))
text(-0.5,0.01,"$Re=100$")
text(-0.5,0.006,"$Re=200$")
text(-0.5,0.001,"$Re=400$")
dev.off()
