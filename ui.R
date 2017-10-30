#==============================================================================
#    ui.R : Genetic Drift Simulator User-Interface
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
library(ggplot2)


#==============================================================================
# shinyUI
#==============================================================================
shinyUI(

	pageWithSidebar(

#------------------------------------------------------------------------------
		headerPanel("Genetic Drift"),
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Sidebar with input
#------------------------------------------------------------------------------
		sidebarPanel(
			wellPanel(
				sliderInput('M',       'Number of gene copies M:', min=10, max=200, value=50, step=5)
				, sliderInput('p0',    'Initial frequency p0 of allele A:', min=0.05, max=0.95, value=0.2, step=0.05)
				, sliderInput('NbGen', 'Number of generations:', min=10, max=1000, value=100, step=10)
				, actionButton('go',   'New Simulation', icon("random"))
			)
			, wellPanel(
				sliderInput('NbRepA', 'Number of repetitions A:', min=1, max=1000, value=1, step=1)
				, sliderInput('NbRepB', 'Number of repetitions B:', min=1, max=100, value=10, step=1)
				, checkboxInput('FixationFlag',  'Show fixation', FALSE)
				, checkboxInput('ExpectedFlag',  'Show expected', FALSE)
			)
		),


#------------------------------------------------------------------------------
# Result Panel
#------------------------------------------------------------------------------
		mainPanel(
			tabsetPanel(
				type="tabs"
				, tabPanel("Plot",                 plotOutput("driftPlotFreq", height = "600px"))
				, tabPanel("Plot Density",         plotOutput("driftPlotFreqDensity", height = "600px"))
				, tabPanel("Mean of P",            plotOutput("driftPlotMeanP", height = "600px"))
				, tabPanel("Variance of P",        plotOutput("driftPlotVarianceP", height = "600px"))
				, tabPanel("Fixation Probability", plotOutput("driftPlotFixationProbability", height = "600px"))
				, tabPanel("Fixation Time",        plotOutput("driftPlotFixationTime", height = "600px"))
			)
		)
#------------------------------------------------------------------------------

	)
)

