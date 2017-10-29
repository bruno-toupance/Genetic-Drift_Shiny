#==============================================================================
#    Genetic-Drift.R : Genetic Drift Simulator
#    Copyright (C) 2017  Bruno Toupance <bruno.toupance@mnhn.fr>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#==============================================================================


require(ggplot2)



#==============================================================================
# SimulateData
#  Simulate Frequency Evolution under Drift for Allele A
#==============================================================================
SimulateData <- function(M=10, p0=0.5, NbGen=15, NbRep=10)
{
#------------------------------------------------------------------------------
	MAT <- matrix(NA, nrow=NbRep, ncol=NbGen+1)
	MAT[, 1] <- p0
#------------------------------------------------------------------------------
	MAT2 <- matrix(0, nrow=M+1, ncol=NbGen+1)
	MAT2[round(p0*M)+1, 1] <- NbRep
#------------------------------------------------------------------------------
	MeanP  <- rep(NA, NbGen+1)
	VarP   <- rep(NA, NbGen+1)
	NbFixA <- rep(NA, NbGen+1)
	NbLosA <- rep(NA, NbGen+1)
	NbPoly <- rep(NA, NbGen+1)
	TimeFA <- c()
	TimeLA <- c()
#------------------------------------------------------------------------------
	jMax <- NbGen+1
	for (Rep in 1:NbRep) {
		i <- Rep
		p <- p0
		Gen <- 1
		while (Gen <= NbGen) {
			j <- Gen+1
			n <- rbinom(1, M, p)
			p <- n/M
			k <- n+1
			if (n == M) {
				TimeFA <- c(TimeFA, Gen)
				MAT[i, j:jMax] <- p
				MAT2[k, j:jMax] <- MAT2[k, j:jMax]+1
				Gen <- NbGen
			} else {
				if (n == 0) {
					TimeLA <- c(TimeLA, Gen)
					MAT[i, j:jMax] <- p
					MAT2[k, j:jMax] <- MAT2[k, j:jMax]+1
					Gen <- NbGen
				} else {
					MAT[i, j] <- p
					MAT2[k, j] <- MAT2[k, j]+1
				}
			}
			Gen <- Gen+1
		}
	}
#------------------------------------------------------------------------------
	NbFixA <- apply( MAT, 2, function(x) { return(sum(x == 1)) } )
	NbLosA <- apply( MAT, 2, function(x) { return(sum(x == 0)) } )
	NbPoly <- apply( MAT, 2, function(x) { return(sum(x != 0 & x != 1)) } )
#------------------------------------------------------------------------------
	MeanP <- apply( MAT, 2, function(x) { return(mean(x)) } )
#------------------------------------------------------------------------------
	if (NbRep > 1) {
		VarP <- apply( MAT, 2, function(x) { return(var(x)) } )
	}
#------------------------------------------------------------------------------
	DATA <- list(MAT=MAT, MAT2=MAT2, M=M, p0=p0, NbGen=NbGen, NbRep=NbRep, NbFixA=NbFixA, NbLosA=NbLosA, NbPoly=NbPoly, MeanP=MeanP, VarP=VarP, TimeFA=TimeFA, TimeLA=TimeLA)
#------------------------------------------------------------------------------
	return(DATA)
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# PlotFreq
#==============================================================================
PlotFreq <- function(DATA)
{
#------------------------------------------------------------------------------
	MAT   <- DATA$MAT
	NbGen <- DATA$NbGen
	NbRep <- DATA$NbRep
#------------------------------------------------------------------------------
	ParBak <- par(no.readonly=TRUE)
#------------------------------------------------------------------------------
	layout(matrix(c(1), nrow=1, ncol=1, byrow=TRUE))
#------------------------------------------------------------------------------
	par(mar=c(bottom=0.5+4, left=0.5+4, top=0.5+0, right=0.5+4))
	BoxY <- c(0, 1)
	BoxX <- c(0, NbGen)
	plot(BoxX, BoxY, type="n", main="", xlab="Time (generations)", ylab="Frequency P", bty="n")
#------------------------------------------------------------------------------
	abline(h=0, col="grey")
	abline(h=1, col="grey")
#------------------------------------------------------------------------------
	Color <- rainbow(NbRep)
	ValX <- 0:NbGen
	for (i in 1:NbRep) {
		ValY <- MAT[i, ]
		lines(ValX, ValY, col=Color[i])
	}
#------------------------------------------------------------------------------
	par(ParBak)
#------------------------------------------------------------------------------
}
#==============================================================================




#==============================================================================
# PlotFreqDensity
#==============================================================================
PlotFreqDensity <- function(DATA)
{
#------------------------------------------------------------------------------
	MAT2  <- DATA$MAT2
	M     <- DATA$M
	NbGen <- DATA$NbGen
	NbRep <- DATA$NbRep
#------------------------------------------------------------------------------
	ParBak <- par(no.readonly=TRUE)
#------------------------------------------------------------------------------
	layout(matrix(c(1), nrow=1, ncol=1, byrow=TRUE))
#------------------------------------------------------------------------------
	par(mar=c(bottom=0.5+4, left=0.5+4, top=0.5+0, right=0.5+4))
	BoxY <- c(0-0.5, M+0.5)
	BoxX <- c(0-0.5, NbGen+0.5)
	plot(BoxX, BoxY, type="n", main="", xlab="Time (generations)", ylab="Allele A count", bty="n")
#------------------------------------------------------------------------------
	Step <- 1/(NbRep+1)
	ScaleArray <- (1-(log10(seq(from=Step, to=1, by=Step))/log10(Step))) 
	ColorArray <- rgb(1, 165/255, 0, alpha=ScaleArray) 
	for (n in 0:M) {
		for (t in 1:NbGen) {
			i <- n+1
			j <- t+1
			Val <- MAT2[i, j]
			Color <- ColorArray[Val+1]
			x1 <- t-0.5
			y1 <- n-0.5
			x2 <- t+0.5
			y2 <- n+0.5
			rect(x1, y1, x2, y2, col=Color, border=NA)
		}
	}
#------------------------------------------------------------------------------
	par(ParBak)
#------------------------------------------------------------------------------
}
#==============================================================================




#==============================================================================
# PlotMeanP
#==============================================================================
PlotMeanP <- function(DATA)
{
#------------------------------------------------------------------------------
	p0    <- DATA$p0
	NbGen <- DATA$NbGen
	MeanP <- DATA$MeanP
#------------------------------------------------------------------------------
	ParBak <- par(no.readonly=TRUE)
#------------------------------------------------------------------------------
	layout(matrix(c(1), nrow=1, ncol=1, byrow=TRUE))
#------------------------------------------------------------------------------
	par(mar=c(bottom=0.5+4, left=0.5+4, top=0.5+0, right=0.5+4))
#------------------------------------------------------------------------------
	BoxY <- c(0, 1)
	BoxX <- c(0, NbGen)
#------------------------------------------------------------------------------
	plot(BoxX, BoxY, type="n", main="", xlab="Time (generations)", ylab="Mean of P", bty="n")
#------------------------------------------------------------------------------
	axis(4, at=c(p0), label=c("p0"), las=2, lwd=0)
#------------------------------------------------------------------------------
	abline(h=0, col="grey")
	abline(h=1, col="grey")
#------------------------------------------------------------------------------
	abline(h=p0, lty=2, col="red")
#------------------------------------------------------------------------------
	ValX <- 0:NbGen
	ValY <- MeanP
	lines(ValX, ValY)
#------------------------------------------------------------------------------
	par(ParBak)
}
#==============================================================================




#==============================================================================
# PlotVarianceP
#==============================================================================
PlotVarianceP <- function(DATA)
{
#------------------------------------------------------------------------------
	M     <- DATA$M
	p0    <- DATA$p0
	NbGen <- DATA$NbGen
	NbRep <- DATA$NbRep
	VarP  <- DATA$VarP
#------------------------------------------------------------------------------
	ParBak <- par(no.readonly=TRUE)
#------------------------------------------------------------------------------
	layout(matrix(c(1), nrow=1, ncol=1, byrow=TRUE))
#------------------------------------------------------------------------------
	par(mar=c(bottom=0.5+4, left=0.5+4, top=0.5+0, right=0.5+4))
#------------------------------------------------------------------------------
	MaxY <- p0*(1-p0)
	if (NbRep > 1) {
		MaxY <- max(c(VarP, MaxY))
	}
	MaxY <- ceiling(MaxY*20)/20
	BoxY <- c(0, MaxY)
	BoxX <- c(0, NbGen)
#------------------------------------------------------------------------------
	plot(BoxX, BoxY, type="n", main="", xlab="Time (generations)", ylab="Variance of P", bty="n")
#------------------------------------------------------------------------------
	axis(4, at=c(p0*(1-p0)), label=c("p0(1-p0)"), las=2, lwd=0)
#------------------------------------------------------------------------------
	abline(h=0, col="grey")
	abline(h=p0*(1-p0), col="grey")
#------------------------------------------------------------------------------
	ValX <- 0:NbGen
	ValY <- p0*(1-p0)*(1-(1-1/M)^ValX)
	lines(ValX, ValY, lty=2, col="red")
#------------------------------------------------------------------------------
	if (NbRep > 1) {
		ValX <- 0:NbGen
		ValY <- VarP
		lines(ValX, ValY)
	} else {
		text(mean(BoxX), mean(BoxY), "NOT AVAILABLE")
	}
#------------------------------------------------------------------------------
	par(ParBak)
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# PlotFixationProbability
#==============================================================================
PlotFixationProbability <- function(DATA)
{
#------------------------------------------------------------------------------
	p0     <- DATA$p0
	NbGen  <- DATA$NbGen
	NbRep  <- DATA$NbRep
	NbFixA <- DATA$NbFixA
	NbLosA <- DATA$NbLosA
	NbPoly <- DATA$NbPoly
#------------------------------------------------------------------------------
	ParBak <- par(no.readonly=TRUE)
#------------------------------------------------------------------------------
	layout(matrix(c(1), nrow=1, ncol=1, byrow=TRUE))
#------------------------------------------------------------------------------
	par(mar=c(bottom=0.5+4, left=0.5+4, top=0.5+0, right=0.5+4))
#------------------------------------------------------------------------------
	BoxY <- c(0, 1)
	BoxX <- c(0, NbGen)
#------------------------------------------------------------------------------
	plot(BoxX, BoxY, type="n", main="", xlab="Time (generations)", ylab="Proportion", bty="n")
#------------------------------------------------------------------------------
	axis(4, at=c(p0, 1-p0), label=c("p0", "1-p0"), las=2, lwd=0)
#------------------------------------------------------------------------------
	abline(h=0, col="grey")
	abline(h=1, col="grey")
	abline(h=p0, lty=2, col="red")
	abline(h=1-p0, lty=2, col="blue")
#------------------------------------------------------------------------------
	ValX <- 0:NbGen
	ValY <- NbFixA/NbRep
	lines(ValX, ValY, col="red")
#------------------------------------------------------------------------------
	ValX <- 0:NbGen
	ValY <- NbLosA/NbRep
	lines(ValX, ValY, col="blue")
#------------------------------------------------------------------------------
	ValX <- 0:NbGen
	ValY <- NbPoly/NbRep
	lines(ValX, ValY, col="black")
#------------------------------------------------------------------------------
	par(ParBak)
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# PlotFixationTime
#==============================================================================
PlotFixationTime <- function(DATA)
{
#------------------------------------------------------------------------------
	M      <- DATA$M
	p0     <- DATA$p0
	NbGen  <- DATA$NbGen
	NbRep  <- DATA$NbRep
	TimeFA <- DATA$TimeFA
	TimeLA <- DATA$TimeLA
#------------------------------------------------------------------------------
	MeanTimeFA <- mean(TimeFA)
	MeanTimeLA <- mean(TimeLA)
#------------------------------------------------------------------------------
	ParBak <- par(no.readonly=TRUE)
#------------------------------------------------------------------------------
	layout(matrix(c(1), nrow=1, ncol=1, byrow=TRUE))
#------------------------------------------------------------------------------
	par(mar=c(bottom=0.5+4, left=0.5+4, top=0.5+0, right=0.5+4))
#------------------------------------------------------------------------------
	DATAFRAME <- data.frame(Time=c(TimeFA, TimeLA), Allele=as.factor(c(rep("A", length(TimeFA)), rep("B", length(TimeLA)))))
	PLOT <- ggplot(DATAFRAME, aes(x=Time, colour=Allele, fill=Allele))
	PLOT <- PLOT + geom_histogram(aes(y=..density..), alpha=0.5, position="identity")
	PLOT <- PLOT + geom_density(alpha=0.2)
	PLOT <- PLOT + labs(x="Time (generations)", y="Density")
	# PLOT <- PLOT + geom_vline(xintercept = c(MeanTimeFA, MeanTimeLA))
	print(PLOT)
#------------------------------------------------------------------------------
	par(ParBak)
#------------------------------------------------------------------------------
}
#==============================================================================


