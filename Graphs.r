setwd("C:/Users/Joe/source/repos/Femke/x64/Release/")

# Run all the simulations:
# 0: No pathogen
# 1:5: R3:R7, no fungicide
# 6:10: R3:R7, with fungicide, absolute resistance
# 11:15: R3:R7, with 2 fungicides, absolute resistance
# 16:20: R3:R7, with fungicide, partial resistance
# 21:25: R3:R7, with 2 fungicides, partial resistance
for(i in 0:25) system(paste("Femke.exe",i),show.output.on.console=FALSE)

# Record the year in which severity increases over 1%
effectiveLife = rep(NA,26)
for(i in 0:25){
	# Read in the year
	Y = read.csv(paste("Sim",i,".csv",sep=""))
	effectiveLife[i+1] = min(which(Y$HAD_loss > 5.0))
}

# Figure 3. S with no pathogen and SLI with a pathogen but no control
png("Figure3.png",width=3.5,height=2.5,units="in",res=300,pointsize=11)
par(mar=c(3,3,0.5,0)+0.1)
#for(i in 0:1) system(paste("Kevin.exe",i),show.output.on.console=FALSE)
D<-read.csv("Dynamics0.csv")
plot(D$H~D$Time,type="l",xlab="Time (degree days)",ylab="Area index",col="gray",lty=2,ann=FALSE)
mtext(side=1,line=2,"Time (degree days)")
mtext(side=2,line=2,"Area index")
D<-read.csv("Dynamics1.csv")
lines(D$H~D$Time,lty=1)
lines((D$latentPyc1 + D$latentAsc1)~D$Time,lty=3)
lines(D$infectious1~D$Time,lty=4)
legend("topleft",col=c(1,"gray",1,1),lty=c(1,2,3,4),c("Healthy","Healthy\n(no disease)","Latent","Infectious"),cex=0.8)
dev.off()

# Figure 4. Severity on R3 and R7 no fungicide, and R3 and R7 with fungicide
png("Figure4.png",width=3.5,height=2.5,units="in",res=300,pointsize=11)
par(mar=c(3,3,0.5,0)+0.1)
#for(i in c(1,5,6,10)) system(paste("Kevin.exe",i),show.output.on.console=FALSE)
D<-read.csv("Dynamics1.csv")
plot(D$Severity~D$Time,type="l",xlab="Time (degree days)",ylab="Severity (%)",ylim=c(0,25),ann=FALSE)
mtext(side=1,line=2,"Time (degree days)")
mtext(side=2,line=2,"Severity (%)")
D<-read.csv("Dynamics5.csv")
lines(D$Severity~D$Time,lty=2)
D<-read.csv("Dynamics6.csv")
lines(D$Severity~D$Time,lty=3)
D<-read.csv("Dynamics10.csv")
lines(D$Severity~D$Time,lty=4)
legend("topleft",lty=1:4,c("R3; no fung", "R7; no fung", "R3; with fung", "R7; with fung"),cex=0.8)
dev.off()

# Figure 5. Resistance frequency over time
png("Figure5.png",width=3.5,height=2.5,units="in",res=300,pointsize=11)
par(mar=c(3,3,0.5,0)+0.1)
D<-read.csv("Sim6.csv")
plot(100*D$Freq2~D$Season,type="l",xlim=c(1,15),xlab="Time (Years)",ylab="Resistance frequency (%)",ylim=c(0,100),ann=FALSE)
mtext(side=1,line=2,"Time (degree days)")
mtext(side=2,line=2,"Resistance frequency (%)")
D<-read.csv("Sim7.csv")
lines(100*D$Freq2~D$Season,lty=2)
D<-read.csv("Sim8.csv")
lines(100*D$Freq2~D$Season,lty=3)
D<-read.csv("Sim9.csv")
lines(100*D$Freq2~D$Season,lty=4)
D<-read.csv("Sim10.csv")
lines(100*D$Freq2~D$Season,lty=5)
legend("topleft",lty=1:5,title="Resistance rating:",c("3","4","5","6","7"),cex=0.8)
dev.off()

# Figure 6. HAD loss over time
png("Figure6.png",width=3.5,height=2.5,units="in",res=300,pointsize=11)
par(mar=c(3,3,0.5,0)+0.1)
D<-read.csv("Sim6.csv")
plot(D$HAD_loss~D$Season,type="l",xlim=c(1,15),xlab="Time (Years)",ylab="HAD loss (%)",ylim=c(0,27),ann=FALSE)
mtext(side=1,line=2,"Time (years)")
mtext(side=2,line=2,"HAD loss (%)")
D<-read.csv("Sim7.csv")
lines(D$HAD_loss~D$Season,lty=2)
D<-read.csv("Sim8.csv")
lines(D$HAD_loss~D$Season,lty=3)
D<-read.csv("Sim9.csv")
lines(D$HAD_loss~D$Season,lty=4)
D<-read.csv("Sim10.csv")
lines(D$HAD_loss~D$Season,lty=5)
legend("topleft",lty=1:5,title="Resistance rating:",c("3","4","5","6","7"),cex=0.8)
dev.off()

# Figure 7. Resistance rating ~ Effective life (years)
png("Figure7.png",width=3.5,height=2.5,units="in",res=300,pointsize=11)
par(mar=c(3,3,0.5,0)+0.1)
plot(effectiveLife[7:11]~seq(3,7),ylim=c(0,100),xlab="Resistance rating",type="o",ylab="Effective life (years)",ann=FALSE)
mtext(side=1,line=2,"Resistance rating")
mtext(side=2,line=2,"Effective life (years)")
lines(effectiveLife[12:16]~seq(3,7),type="o",lty=2)
lines(effectiveLife[17:21]~seq(3,7),type="o",lty=3)
lines(effectiveLife[22:26]~seq(3,7),type="o",lty=4)
legend("topright",lty=c(1,2,3,4),c("AR; Solo","AR; Mix", "PR; Solo", "PR; Mix"),cex=0.8)
dev.off()
