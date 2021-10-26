#==============================================================================
#    ui.R : Genetic Drift Simulator User-Interface
#    Copyright (C) 2021  Bruno Toupance <bruno.toupance@mnhn.fr>
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

#----------------------- Disable data export using 'shinyFiles'
# library(shinyFiles)
#----------------------- Disable data export using 'shinyFiles'


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
				  checkboxInput(inputId = 'DipFlag', label = 'Diploid', TRUE)
				, numericInput(inputId = 'N', label = 'Population size - integer:', value = 50)
				# , numericInput(inputId = 'N', label = 'Population size - integer:', value = 16)
				, numericInput(inputId = 'NbGen', label = 'Number of generations - integer:', value = 100)
				# , numericInput(inputId = 'NbGen', label = 'Number of generations - integer:', value = 20)
				, numericInput(inputId = 'p0', label = 'Initial frequency of allele A - numeric [0, 1]:', value = 0.2)
				# , numericInput(inputId = 'p0', label = 'Initial frequency p(0) of allele A - numeric [0, 1]:', value = 0.5, min = 0, max = 1, step = 0.1)
				, numericInput(inputId = 'NbRep', label = 'Number of repetitions - integer:', value = 1)
				# , numericInput(inputId = 'NbRep', label = 'Number of repetitions - integer:', value = 105)
				, actionButton(inputId = 'go', label = 'New Simulation', icon("random"))
				, actionButton(inputId = 'mul2', label = 'x 2')
				, actionButton(inputId = 'div2', label = 'x 1/2')
			)
		),


#------------------------------------------------------------------------------
# Result Panel
#------------------------------------------------------------------------------
		mainPanel(
			tabsetPanel(type = "tabs"
				, tabPanel(title = "Plot"
					, checkboxInput(inputId = 'FixFlag', label = 'Show fixation', FALSE)
					, plotOutput(outputId = "driftPlotFreq", height = "600px")
					)
				, tabPanel(title = "Plot Density"
					, splitLayout(
						checkboxInput(inputId = 'CntFlag', label = 'Show count', FALSE)
						# , actionButton(inputId = 'export_count', label = 'Export')
						# , shinyFilesButton('files', label = 'Export', title = 'Please select a file', multiple = FALSE)
#----------------------- Disable data export using 'shinyFiles'
#						, shinySaveButton(id = "CntSaveBtn", label = "Save file", title = "Save file as ...", filetype = list(txt = "txt"))
#----------------------- Disable data export using 'shinyFiles'
						)
					, plotOutput(outputId = "driftPlotFreqDensity", height = "600px")
					)
					
				, tabPanel(title = "Mean of P"
					, checkboxInput(inputId = 'ExpMeanFlag', label = 'Show expected', FALSE)
					, plotOutput(outputId = "driftPlotMeanP", height = "600px")
					)
				, tabPanel(title = "Variance of P"
					, checkboxInput(inputId = 'ExpVarFlag', label = 'Show expected', FALSE)
					, plotOutput(outputId = "driftPlotVarianceP", height = "600px")
					)
				, tabPanel(title = "Fixation Probability"
					, checkboxInput(inputId = 'ExpFixFlag', label = 'Show expected', FALSE)
					, plotOutput(outputId = "driftPlotFixationProbability", height = "600px")
					)
				, tabPanel(title = "Fixation Time"
					, checkboxInput(inputId = 'ExpTimeFlag', label = 'Show expected', FALSE)
					, plotOutput(outputId = "driftPlotFixationTime", height = "600px")
					)
			)
		)
#------------------------------------------------------------------------------

	)
)

