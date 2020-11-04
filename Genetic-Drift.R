#==============================================================================
#    Genetic-Drift.R : Genetic Drift Simulator
#    Copyright (C) 2020  Bruno Toupance <bruno.toupance@mnhn.fr>
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
require(colorspace)
# require(scales)


#==============================================================================
# check_parameters
#==============================================================================
check_parameters <- function(N=10, p0=0.5, NbGen=15, NbRep=10) 
{
#------------------------------------------------------------------------------
	Flag <- TRUE
	Msg <- ""
#------------------------------------------------------------------------------
	if (is.numeric(N)) {
		if (N == round(N)) {
			if ( N < 1 ) {
				Flag <- FALSE
				Msg <- sprintf("%s\n%s", Msg, "FAIL: [ N ] out of bounds")
			}
		} else {
			Flag <- FALSE
			Msg <- sprintf("%s\n%s", Msg, "FAIL: [ N ] not integer")
		}
	} else {
		Flag <- FALSE
		Msg <- sprintf("%s\n%s", Msg, "FAIL: [ N ] not numeric")
	}
#------------------------------------------------------------------------------
	if (is.numeric(p0)) {
		if ( (p0 < 0) | (p0 > 1) ) {
			Flag <- FALSE
			Msg <- sprintf("%s\n%s", Msg, "FAIL: [ p(0) ] out of bounds")
		}
	} else {
		Flag <- FALSE
		Msg <- sprintf("%s\n%s", Msg, "FAIL: [ p(0) ] not numeric")
	}
#------------------------------------------------------------------------------
	if (is.numeric(NbGen)) {
		if (NbGen == round(NbGen)) {
			if ( NbGen < 1 ) {
				Flag <- FALSE
				Msg <- sprintf("%s\n%s", Msg, "FAIL: [ NbGen ] out of bounds")
			}
		} else {
			Flag <- FALSE
			Msg <- sprintf("%s\n%s", Msg, "FAIL: [ NbGen ] not integer")
		}
	} else {
		Flag <- FALSE
		Msg <- sprintf("%s\n%s", Msg, "FAIL: [ NbGen ] not numeric")
	}
#------------------------------------------------------------------------------
	if (is.numeric(NbRep)) {
		if (NbRep == round(NbRep)) {
			if ( NbRep < 1 ) {
				Flag <- FALSE
				Msg <- sprintf("%s\n%s", Msg, "FAIL: [ NbRep ] out of bounds")
			}
		} else {
			Flag <- FALSE
			Msg <- sprintf("%s\n%s", Msg, "FAIL: [ NbRep ] not integer")
		}
	} else {
		Flag <- FALSE
		Msg <- sprintf("%s\n%s", Msg, "FAIL: [ NbRep ] not numeric")
	}
#------------------------------------------------------------------------------
	return(list(Msg=Msg, Flag=Flag, N=N, p0=p0, NbGen=NbGen, NbRep=NbRep))
#------------------------------------------------------------------------------
}





#==============================================================================
# SimulateData
#  Simulate Frequency Evolution under Drift for Allele A
#==============================================================================
SimulateData <- function(N=10, p0=0.5, NbGen=15, NbRep=10, DipFlag=TRUE) 
{
#------------------------------------------------------------------------------
	ParamChecking <- check_parameters(N, p0, NbGen, NbRep)
	# print(ParamChecking)
#------------------------------------------------------------------------------
	if (ParamChecking$Flag) {
		if (DipFlag) {
			NC <- 2*N
		} else {
			NC <- N
		}
		FreqMat <- matrix(NA, nrow=NbRep, ncol=NbGen+1)
		FreqMat[, 1] <- p0

		CountMat <- matrix(0, nrow=NC+1, ncol=NbGen+1)
		CountMat[round(p0*NC)+1, 1] <- NbRep

		MeanP  <- rep(NA, NbGen+1)
		VarP   <- rep(NA, NbGen+1)
		NbFixA <- rep(NA, NbGen+1)
		NbLosA <- rep(NA, NbGen+1)
		NbPoly <- rep(NA, NbGen+1)
		TimeFA <- c()
		TimeLA <- c()

		jMax <- NbGen+1
		for (Rep in 1:NbRep) {
			i <- Rep
			p <- p0
			Gen <- 1
			while (Gen <= NbGen) {
				j <- Gen+1
				n <- rbinom(1, NC, p)
				p <- n/NC
				k <- n+1
				if (n == NC) {
					TimeFA <- c(TimeFA, Gen)
					FreqMat[i, j:jMax] <- p
					CountMat[k, j:jMax] <- CountMat[k, j:jMax]+1
					Gen <- NbGen
				} else {
					if (n == 0) {
						TimeLA <- c(TimeLA, Gen)
						FreqMat[i, j:jMax] <- p
						CountMat[k, j:jMax] <- CountMat[k, j:jMax]+1
						Gen <- NbGen
					} else {
						FreqMat[i, j] <- p
						CountMat[k, j] <- CountMat[k, j]+1
					}
				}
				Gen <- Gen+1
			}
		}

		NbFixA <- apply(FreqMat, 2, function(x) { return(sum(x == 1)) } )
		NbLosA <- apply(FreqMat, 2, function(x) { return(sum(x == 0)) } )
		NbPoly <- apply(FreqMat, 2, function(x) { return(sum(x != 0 & x != 1)) } )

		MeanP <- apply(FreqMat, 2, function(x) { return(mean(x)) } )
		if (NbRep > 1) {
			VarP <- apply(FreqMat, 2, function(x) { return(var(x)) } )
		}

		SimData <- list(FreqMat=FreqMat, CountMat=CountMat, NC=NC, p0=p0, 
			NbGen=NbGen, NbRep=NbRep, 
			NbFixA=NbFixA, NbLosA=NbLosA, NbPoly=NbPoly, 
			MeanP=MeanP, VarP=VarP, 
			TimeFA=TimeFA, TimeLA=TimeLA, Param=ParamChecking)
#------------------------------------------------------------------------------
	} else {
		SimData <- list(Param=ParamChecking)
	}
#------------------------------------------------------------------------------
	return(SimData)
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# PlotError
#==============================================================================
PlotError <- function(Msg)
{
#------------------------------------------------------------------------------
	plot(c(0, 1), c(0, 1), type="n", xlab="", ylab="", main="", 
		xaxt="n", yaxt="n", bty="n")
	ErrMsg <- sprintf("ERROR: check parameter values%s", Msg)
	text(0.5, 0.5, ErrMsg, col="red", adj=0.5)
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# PlotFreq
#==============================================================================
PlotFreq <- function(SimData, FixFlag=FALSE)
{
#------------------------------------------------------------------------------
	if (SimData$Param$Flag) {
		FreqMat <- SimData$FreqMat
		NbGen <- SimData$NbGen
		NbRep <- SimData$NbRep

		ParBak <- par(no.readonly=TRUE)

		layout(matrix(c(1), nrow=1, ncol=1, byrow=TRUE))

		par(mar=c(bottom=0.5+4, left=0.5+4, top=0.5+0, right=0.5+5))
		BoxY <- c(0, 1)
		BoxX <- c(0, NbGen)
		plot(BoxX, BoxY, type="n", main="", xlab="Time (generations)", 
			ylab="Frequency P", bty="n")

		abline(h=0, col="gray")
		abline(h=1, col="gray")
	
		# Color <- hue_pal()(NbRep)  ## 'scales' package
		# Color <- rainbow(NbRep)
		Color <- rainbow_hcl(NbRep)   ## 'colorspace' package
		Color_AB <- rainbow_hcl(2)   ## 'colorspace' package
		Color_FA <- Color_AB[1]
		Color_LA <- Color_AB[2]
		
		ValX <- 0:NbGen
		FixA <- FreqMat[, NbGen+1] == 1
		FixB <- FreqMat[, NbGen+1] == 0
		Poly <- ! (FixA | FixB) 
		for (i in 1:NbRep) {
			ValY <- FreqMat[i, ]
			if (FixFlag) {
				Pos <- which((ValY != 0) & (ValY != 1))
				Pos0 <- which(ValY == 0)
				Pos1 <- which(ValY == 1)
				if (length(Pos0) > 1) {
					Pos0 <- Pos0[1]
					Pos <- c(Pos, Pos0)
					lines(ValX[Pos], ValY[Pos], col=Color_LA)
					points(ValX[Pos0], ValY[Pos0], pch=20, col=Color_LA)
				} else {
					if (length(Pos1) > 1) {
						Pos1 <- Pos1[1]
						Pos <- c(Pos, Pos1)
						lines(ValX[Pos], ValY[Pos], col=Color_FA)
						points(ValX[Pos1], ValY[Pos1], pch=20, col=Color_FA)
					} else {
						lines(ValX, ValY, col="black")
					}
				}
			} else {
				lines(ValX, ValY, col=Color[i])
			}
		}
		ValX <- 0
		ValY <- SimData$p0
		points(ValX, ValY, pch=15, col="black", cex=1)

		if (FixFlag) {
			FixB_Label <- sprintf("Loss=%d", sum(FixB))
			Poly_Label <- sprintf("Polym=%d", sum(Poly))
			FixA_Label <- sprintf("Fixed=%d", sum(FixA))
			axis(4, at=c(0.0), label=FixB_Label, las=2, lwd=0, col.axis=Color_LA, 
				cex.axis=0.7)
			axis(4, at=c(0.5), label=Poly_Label, las=2, lwd=0, col.axis="black", 
				cex.axis=0.7)
			axis(4, at=c(1.0), label=FixA_Label, las=2, lwd=0, col.axis=Color_FA, 
				cex.axis=0.7)
			if (sum(SimData$TimeFA) > 0) {
				rug(SimData$TimeFA, side=3, col=Color_FA)
				rug(mean(SimData$TimeFA), side=3, col="black", lwd=2)
			}
			if (sum(SimData$TimeLA) > 0) {
				rug(SimData$TimeLA, side=1, col=Color_LA)
				rug(mean(SimData$TimeLA), side=1, col="black", lwd=2)
			}
		}

		par(ParBak)
#------------------------------------------------------------------------------
	} else {
		PlotError(SimData$Param$Msg)
	}
#------------------------------------------------------------------------------
}
#==============================================================================




#==============================================================================
# PlotFreqDensity
#==============================================================================
PlotFreqDensity <- function(SimData, CntFlag=FALSE) 
{
#------------------------------------------------------------------------------
	if (SimData$Param$Flag) {
		CountMat <- SimData$CountMat
		NC       <- SimData$NC
		NbGen    <- SimData$NbGen
		NbRep    <- SimData$NbRep

		ParBak <- par(no.readonly=TRUE)

		layout(matrix(c(1), nrow=1, ncol=1, byrow=TRUE))

		par(mar=c(bottom=0.5+4, left=0.5+4, top=0.5+0, right=0.5+5))
		BoxY <- c(0-0.5, NC+0.5)
		BoxX <- c(0-0.5, NbGen+0.5)
		plot(BoxX, BoxY, type="n", main="", xlab="Time (generations)", 
			ylab="Allele A count", bty="n")

		Step <- 1/(NbRep+1)
		ScaleArray <- (1-(log10(seq(from=Step, to=1, by=Step))/log10(Step))) 
		ColorArray <- rgb(1, 165/255, 0, alpha=ScaleArray) 
		for (n in 0:NC) {
			for (t in 1:NbGen) {
				i <- n+1
				j <- t+1
				Val <- CountMat[i, j]
				Color <- ColorArray[Val+1]
				x1 <- t-0.5
				y1 <- n-0.5
				x2 <- t+0.5
				y2 <- n+0.5
				rect(x1, y1, x2, y2, col=Color, border=NA)
				if (CntFlag) {
					text(t, n, Val)
				}
			}
		}

		par(ParBak)
#------------------------------------------------------------------------------
	} else {
		PlotError(SimData$Param$Msg)
	}
#------------------------------------------------------------------------------
}
#==============================================================================




#==============================================================================
# PlotMeanP
#==============================================================================
PlotMeanP <- function(SimData, ExpFlag=FALSE) 
{
#------------------------------------------------------------------------------
	if (SimData$Param$Flag) {
		p0    <- SimData$p0
		NbGen <- SimData$NbGen
		MeanP <- SimData$MeanP

		ParBak <- par(no.readonly=TRUE)

		layout(matrix(c(1), nrow=1, ncol=1, byrow=TRUE))

		par(mar=c(bottom=0.5+4, left=0.5+4, top=0.5+0, right=0.5+5))

		BoxY <- c(0, 1)
		BoxX <- c(0, NbGen)

		plot(BoxX, BoxY, type="n", main="", xlab="Time (generations)", 
			ylab="Mean of P", bty="n")

		abline(h=seq(from=0, to=1, by=0.1), col="gray", lty=2)

		if (ExpFlag) {
			axis(4, at=c(p0), label=c("p(0)"), las=2, lwd=0, cex.axis=0.7)
			# axis(4, at=p0, label=expression(p[0]), las=2, lwd=0)
			abline(h=p0, lty=2, col="red")
		}

		abline(h=0, col="gray")
		abline(h=1, col="gray")

		ValX <- 0:NbGen
		ValY <- MeanP
		lines(ValX, ValY)

		if (NbGen < 50) {
			points(ValX, ValY, pch=15, col="gray", cex=1)
		}
		
		par(ParBak)
#------------------------------------------------------------------------------
	} else {
		PlotError(SimData$Param$Msg)
	}
#------------------------------------------------------------------------------
}
#==============================================================================




#==============================================================================
# PlotVarianceP
#==============================================================================
PlotVarianceP <- function(SimData, ExpFlag=FALSE) 
{
#------------------------------------------------------------------------------
	if (SimData$Param$Flag) {
		NC     <- SimData$NC
		p0    <- SimData$p0
		NbGen <- SimData$NbGen
		NbRep <- SimData$NbRep
		VarP  <- SimData$VarP

		ParBak <- par(no.readonly=TRUE)

		layout(matrix(c(1), nrow=1, ncol=1, byrow=TRUE))

		par(mar=c(bottom=0.5+4, left=0.5+4, top=0.5+0, right=0.5+5))

		MaxY <- p0*(1-p0)
		if (NbRep > 1) {
			MaxY <- max(c(VarP, MaxY))
		}
		MaxY <- ceiling(MaxY*20)/20
		BoxY <- c(0, MaxY)
		BoxX <- c(0, NbGen)

		plot(BoxX, BoxY, type="n", main="", xlab="Time (generations)", 
			ylab="Variance of P", bty="n")

		if (ExpFlag) {
			axis(4, at=c(p0*(1-p0)), label=c("p(0)(1-p(0))"), las=2, lwd=0, 
				cex.axis=0.7)
			# axis(4, at=p0*(1-p0), label=expression(p[0](1-p[0])), las=2, lwd=0)
			abline(h=p0*(1-p0), col="gray")
			ValX <- 0:NbGen
			ValY <- p0*(1-p0)*(1-(1-1/NC)^ValX)
			lines(ValX, ValY, lty=2, col="red")
		}

		abline(h=0, col="gray")

		if (NbRep > 1) {
			ValX <- 0:NbGen
			ValY <- VarP
			lines(ValX, ValY)

			if (NbGen < 50) {
				points(ValX, ValY, pch=15, col="gray", cex=1)
			}
		
		} else {
			text(mean(BoxX), mean(BoxY), "NOT AVAILABLE", col="red")
		}

		par(ParBak)
#------------------------------------------------------------------------------
	} else {
		PlotError(SimData$Param$Msg)
	}
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# PlotFixationProbability
#==============================================================================
PlotFixationProbability <- function(SimData, ExpFlag=FALSE) 
{
#------------------------------------------------------------------------------
	if (SimData$Param$Flag) {
		p0     <- SimData$p0
		NbGen  <- SimData$NbGen
		NbRep  <- SimData$NbRep
		NbFixA <- SimData$NbFixA
		NbLosA <- SimData$NbLosA
		NbPoly <- SimData$NbPoly

		ParBak <- par(no.readonly=TRUE)

		Color_AB <- rainbow_hcl(2)  ## 'colorspace' package
		Color_FA <- Color_AB[1]
		Color_LA <- Color_AB[2]

		layout(matrix(c(1), nrow=1, ncol=1, byrow=TRUE))

		par(mar=c(bottom=0.5+4, left=0.5+4, top=0.5+0, right=0.5+5))

		BoxY <- c(0, 1)
		BoxX <- c(0, NbGen)

		plot(BoxX, BoxY, type="n", main="", xlab="Time (generations)", 
			ylab="Proportion", bty="n")

		abline(h=seq(from=0, to=1, by=0.1), col="gray", lty=2)
		
		if (ExpFlag) {
			abline(h=p0, lty=2, col=Color_FA)
			abline(h=1-p0, lty=2, col=Color_LA)
			if (p0 == 0.5) {
				axis(4, at=p0, label="p(0)=1-p(0)", las=2, lwd=0, cex.axis=0.7)
			} else {
				axis(4, at=p0, label="p(0)", las=2, lwd=0, cex.axis=0.7)
				axis(4, at=1-p0, label="1-p(0)", las=2, lwd=0, cex.axis=0.7)
			}
			# axis(4, at=c(p0, 1-p0), label=c(expression(p[0]), 
				# expression(1-p[0])), las=2, lwd=0)
		}

		abline(h=0, col="gray")
		abline(h=1, col="gray")

		ValX <- 0:NbGen
		ValY <- NbFixA/NbRep
		lines(ValX, ValY, col=Color_FA)

		ValX <- 0:NbGen
		ValY <- NbLosA/NbRep
		lines(ValX, ValY, col=Color_LA)

		ValX <- 0:NbGen
		ValY <- NbPoly/NbRep
		lines(ValX, ValY, col="black")
		
		legend("topright", legend=c("Polym", "A fixed", "a fixed"), 
			col=c("black", Color_FA, Color_LA), lty=1, bg="white")

		par(ParBak)
#------------------------------------------------------------------------------
	} else {
		PlotError(SimData$Param$Msg)
	}
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# PlotFixationTime
#==============================================================================
PlotFixationTime <- function(SimData, ExpFlag=FALSE)
{
#------------------------------------------------------------------------------
	if (SimData$Param$Flag) {
		NC     <- SimData$NC
		p0     <- SimData$p0
		NbGen  <- SimData$NbGen
		NbRep  <- SimData$NbRep
		TimeFA <- SimData$TimeFA
		TimeLA <- SimData$TimeLA

		MeanTimeFA <- mean(TimeFA)
		MeanTimeLA <- mean(TimeLA)

		# Color_AB <- hue_pal()(2)  ## 'scales' package
		Color_AB <- rainbow_hcl(2)  ## 'colorspace' package

		ParBak <- par(no.readonly=TRUE)

		layout(matrix(c(1), nrow=1, ncol=1, byrow=TRUE))

		par(mar=c(bottom=0.5+4, left=0.5+4, top=0.5+0, right=0.5+5))

		ExpTimeFA <- -2*NC*(1-p0)/p0*log(1-p0)
		ExpTimeLA <- -2*NC*(p0)/(1-p0)*log(p0)


		A_label <- sprintf("A\nObserved mean = %5.1f", MeanTimeFA)
		B_label <- sprintf("B\nObserved mean = %5.1f", MeanTimeLA)
		if (ExpFlag) {
			A_label <- sprintf("%s\nExpected mean = %5.1f", A_label, ExpTimeFA)
			B_label <- sprintf("%s\nExpected mean = %5.1f", B_label, ExpTimeLA)
		}
	
		if (length(TimeFA) > 0 | length(TimeLA) > 0) {
			PlotData <- data.frame(Time=c(TimeFA, TimeLA), 
				Allele=as.factor(c(rep(A_label, length(TimeFA)), 
				rep(B_label, length(TimeLA)))))
			Plot <- ggplot(PlotData, aes(x=Time, colour=Allele, fill=Allele))
			Plot <- Plot + geom_histogram(aes(y=..density..), alpha=0.5, 
				position="identity")
			Plot <- Plot + xlim(c(0, NbGen))
			Plot <- Plot + geom_density(alpha=0.2)
			Plot <- Plot + labs(x="Time (generations)", y="Density")
			Plot <- Plot + geom_vline(xintercept=c(MeanTimeFA, MeanTimeLA), 
				colour=Color_AB, size=1)
			if (ExpFlag) {
				Plot <- Plot + geom_vline(xintercept=c(ExpTimeFA, ExpTimeLA), 
					colour=Color_AB, linetype="dotted", size=1)
			}
			Plot <- Plot + theme(legend.justification=c(1, 1), 
				legend.position=c(1, 1))

			print(Plot)
		} else {
			plot(c(0, 1), c(0, 1), type="n", xlab="", ylab="", main="", 
				xaxt="n", yaxt="n", bty="n")
			ErrMsg <- "NO FIXATION"
			text(0.5, 0.5, ErrMsg, col="red", adj=0.5)
		}

		par(ParBak)
#------------------------------------------------------------------------------
	} else {
		PlotError(SimData$Param$Msg)
	}
#------------------------------------------------------------------------------
}
#==============================================================================


