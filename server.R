#==============================================================================
#    server.R : Genetic Drift Simulator Server
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


library(shiny)

source("Genetic-Drift.R")

#==============================================================================
# shinyServer
#==============================================================================
shinyServer(
	function(input, output) {

#------------------------------------------------------------------------------
		DATA <- eventReactive(input$go, {
			SimulateData(input$M, input$p0, input$NbGen, input$NbRep)
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$driftPlotFreq <- renderPlot({
			PLOT <- PlotFreq(DATA())
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$driftPlotFreqDensity <- renderPlot({
			PLOT <- PlotFreqDensity(DATA())
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$driftPlotMeanP <- renderPlot({
			PLOT <- PlotMeanP(DATA())
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$driftPlotVarianceP <- renderPlot({
			PLOT <- PlotVarianceP(DATA())
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$driftPlotFixationProbability <- renderPlot({
			PLOT <- PlotFixationProbability(DATA())
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$driftPlotFixationTime <- renderPlot({
			PLOT <- PlotFixationTime(DATA())
		})
#------------------------------------------------------------------------------


	}
)
#==============================================================================
