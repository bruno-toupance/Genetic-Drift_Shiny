#==============================================================================
#    server.R : Genetic Drift Simulator Server
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


library(shiny)

source("Genetic-Drift.R")

#==============================================================================
# shinyServer
#==============================================================================
shinyServer(
	function(input, output, session) {

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
		observeEvent(input$mul2, {
			updateNumericInput(session, "N", value=input$N*2)
			updateNumericInput(session, "NbGen", value=input$NbGen*2)
		})
#------------------------------------------------------------------------------
		observeEvent(input$div2, {
			updateNumericInput(session, "N", value=input$N %/% 2)
			updateNumericInput(session, "NbGen", value=input$NbGen %/% 2)
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$driftPlotFreq <- renderPlot({
			Plot <- PlotFreq(DriftData(), FixFlag=input$FixFlag)
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$driftPlotFreqDensity <- renderPlot({
			Plot <- PlotFreqDensity(DriftData(), CntFlag=input$CntFlag)
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$driftPlotMeanP <- renderPlot({
			Plot <- PlotMeanP(DriftData(), ExpFlag=input$ExpFlag)
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$driftPlotVarianceP <- renderPlot({
			Plot <- PlotVarianceP(DriftData(), ExpFlag=input$ExpFlag)
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$driftPlotFixationProbability <- renderPlot({
			Plot <- PlotFixationProbability(DriftData(), ExpFlag=input$ExpFlag)
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$driftPlotFixationTime <- renderPlot({
			Plot <- PlotFixationTime(DriftData(), ExpFlag=input$ExpFlag)
		})
#------------------------------------------------------------------------------


	}
)
#==============================================================================
