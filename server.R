#==============================================================================
#    server.R : Genetic Drift Simulator Server
#    Copyright (C) 2019  Bruno Toupance <bruno.toupance@mnhn.fr>
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


library(shiny)

source("Genetic-Drift.R")

#==============================================================================
# shinyServer
#==============================================================================
shinyServer(
	function(input, output) {

#------------------------------------------------------------------------------
		DriftData <- reactive({
			Tmp <- input$go
			return(SimulateData(
				N=input$N, 
				p0=input$p0, 
				NbGen=input$NbGen, 
				NbRep=input$NbRep,
				DipFlag=input$DipFlag))
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$driftPlotFreq <- renderPlot({
			Plot <- PlotFreq(DriftData(), input$FixFlag)
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$driftPlotFreqDensity <- renderPlot({
			Plot <- PlotFreqDensity(DriftData())
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$driftPlotMeanP <- renderPlot({
			Plot <- PlotMeanP(DriftData(), input$ExpFlag)
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$driftPlotVarianceP <- renderPlot({
			Plot <- PlotVarianceP(DriftData(), input$ExpFlag)
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$driftPlotFixationProbability <- renderPlot({
			Plot <- PlotFixationProbability(DriftData(), input$ExpFlag)
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$driftPlotFixationTime <- renderPlot({
			PLOT <- PlotFixationTime(DriftData())
		})
#------------------------------------------------------------------------------


	}
)
#==============================================================================
