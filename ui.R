#==============================================================================
#    ui.R : Genetic Drift Simulator User-Interface
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
				  checkboxInput(inputId='DipFlag',      label='Diploid', TRUE)
				, numericInput(inputId='N',             label='Population size - integer:', value=50)
				# , numericInput(inputId='N',             label='Population size - integer:', value=16)
				, numericInput(inputId='NbGen',         label='Number of generations - integer:', value=100)
				# , numericInput(inputId='NbGen',         label='Number of generations - integer:', value=20)
				, numericInput(inputId='p0',            label='Initial frequency p(0) of allele A - numeric [0, 1]:', value=0.2)
				# , numericInput(inputId='p0',            label='Initial frequency p(0) of allele A - numeric [0, 1]:', value=0.5, min=0, max=1, step=0.1)
				, numericInput(inputId='NbRep',         label='Number of repetitions - integer:', value=1)
				# , numericInput(inputId='NbRep',         label='Number of repetitions - integer:', value=105)
				, actionButton(inputId='go',            label='New Simulation', icon("random"))
				, actionButton(inputId='mul2',            label='x 2')
				, actionButton(inputId='div2',            label='x 1/2')
			)
			, wellPanel(
				  checkboxInput(inputId='FixFlag',  label='Show fixation', FALSE)
				, checkboxInput(inputId='ExpFlag',  label='Show expected', FALSE)
				, checkboxInput(inputId='CntFlag',  label='Show count', FALSE)
			)
		),


#------------------------------------------------------------------------------
# Result Panel
#------------------------------------------------------------------------------
		mainPanel(
			tabsetPanel(
				  type="tabs"
				, tabPanel(title="Plot",                 plotOutput(outputId="driftPlotFreq",                height = "600px"))
				, tabPanel(title="Plot Density",         plotOutput(outputId="driftPlotFreqDensity",         height = "600px"))
				, tabPanel(title="Mean of P",            plotOutput(outputId="driftPlotMeanP",               height = "600px"))
				, tabPanel(title="Variance of P",        plotOutput(outputId="driftPlotVarianceP",           height = "600px"))
				, tabPanel(title="Fixation Probability", plotOutput(outputId="driftPlotFixationProbability", height = "600px"))
				, tabPanel(title="Fixation Time",        plotOutput(outputId="driftPlotFixationTime",        height = "600px"))
			)
		)
#------------------------------------------------------------------------------

	)
)

