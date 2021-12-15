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


library("shiny")


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
				checkboxInput(inputId = 'diploid_flag', label = 'Diploid', value = TRUE), 
				numericInput(inputId = 'nb_ind', label = 'Population size:', min = 1, value = 50), 
				numericInput(inputId = 'ini_p', label = 'Initial frequency of allele A:', min = 0, max = 1, step = 0.1, value = 0.2), 
				numericInput(inputId = 'nb_gen', label = 'Number of generations:', min = 1, value = 100), 
				numericInput(inputId = 'nb_rep', label = 'Number of repetitions:', min = 1, value = 1), 
				actionButton(inputId = 'go', label = 'New Simulation', icon("random")), 
				actionButton(inputId = 'mul2', label = 'x 2'), 
				actionButton(inputId = 'div2', label = 'x 1/2')
			)
		),


#------------------------------------------------------------------------------
# Result Panel
#------------------------------------------------------------------------------
		mainPanel(
			tabsetPanel(
				type = "tabs",
				
				tabPanel(
					title = "Plot", 
					checkboxInput(inputId = 'fix_flag', label = 'Show fixation', value = FALSE), 
					plotOutput(outputId = "driftPlotFreq", height = "600px")
				),

				tabPanel(
					title = "Plot Density",
					splitLayout(
						checkboxInput(inputId = 'count_flag', label = 'Show count', value = FALSE)
					), 
					plotOutput(outputId = "driftPlotFreqDensity", height = "600px")
				),
					
				tabPanel(
					title = "Count Table",
					splitLayout(
						downloadButton(outputId = 'export_count', label = 'Download')
					), 
					tableOutput(outputId = "count_table")
				),
					
				tabPanel(
					title = "Mean of P", 
					checkboxInput(inputId = 'expected_mean_flag', label = 'Show expected', value = FALSE), 
					plotOutput(outputId = "driftPlotMeanP", height = "600px")
				),

				tabPanel(
					title = "Variance of P", 
					checkboxInput(inputId = 'expected_var_flag', label = 'Show expected', value = FALSE), 
					plotOutput(outputId = "driftPlotVarianceP", height = "600px")
				),

				tabPanel(
					title = "Fixation Probability", 
					checkboxInput(inputId = 'expected_prob_flag', label = 'Show expected', value = FALSE), 
					plotOutput(outputId = "driftPlotFixationProbability", height = "600px")
				),

				tabPanel(
					title = "Fixation Time", 
					checkboxInput(inputId = 'expected_time_flag', label = 'Show expected', value = FALSE),
					plotOutput(outputId = "driftPlotFixationTime", height = "600px")
				)

			)
		)
#------------------------------------------------------------------------------

	)
)

