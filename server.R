#==============================================================================
#    server.R : Genetic Drift Simulator Server
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

source("Genetic-Drift.R")


#==============================================================================
# shinyServer
#==============================================================================
shinyServer(
	function(input, output, session) {

#------------------------------------------------------------------------------
		drift_data <- reactive(
			{
				tmp <- input$go
				return(
					SimulateData(
						nb_ind = input$nb_ind, 
						ini_p = input$ini_p , 
						nb_gen = input$nb_gen, 
						nb_rep = input$nb_rep, 
						diploid_flag = input$diploid_flag
					)
				)
			}
		)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		observeEvent(input$mul2, {
			updateNumericInput(session, "nb_ind", value = input$nb_ind * 2)
			updateNumericInput(session, "nb_gen", value = input$nb_gen * 2)
		})
#------------------------------------------------------------------------------
		observeEvent(input$div2, {
			updateNumericInput(session, "nb_ind", value = input$nb_ind %/% 2)
			updateNumericInput(session, "nb_gen", value = input$nb_gen %/% 2)
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$export_count <- downloadHandler(
				filename = function() {
					alpha_num <- c(0:9, LETTERS[1:6])
					rnd_num <- paste(sample(alpha_num, size = 8, replace = TRUE), collapse = "")
					file_path <- sprintf("drift_data_%s_%s.txt", Sys.Date(), rnd_num)
					return(file_path)
				},
				content = function(file_path) {
					write.table(get_count_df(drift_data()), 
						file = file_path, sep = "\t", 
						quote = FALSE, row.names = FALSE)
				}
		)
#------------------------------------------------------------------------------
		output$count_table <- renderTable(
			{
				get_count_df(drift_data())
			}, 
			digits = 0
		)
#------------------------------------------------------------------------------
		output$driftPlotFreq <- renderPlot(
			{
				Plot <- PlotFreq(drift_data(), fix_flag = input$fix_flag)
			}
		)
#------------------------------------------------------------------------------
		output$driftPlotFreqDensity <- renderPlot(
			{
				Plot <- PlotFreqDensity(drift_data(), count_flag = input$count_flag)
			}
		)
#------------------------------------------------------------------------------
		output$driftPlotMeanP <- renderPlot(
			{
				Plot <- PlotMeanP(drift_data(), expected_mean_flag = input$expected_mean_flag)
			}
		)
#------------------------------------------------------------------------------
		output$driftPlotVarianceP <- renderPlot(
			{
				Plot <- PlotVarianceP(drift_data(), expected_var_flag = input$expected_var_flag)
			}
		)
#------------------------------------------------------------------------------
		output$driftPlotFixationProbability <- renderPlot(
			{
				Plot <- PlotFixationProbability(drift_data(), expected_prob_flag = input$expected_prob_flag)
			}
		)
#------------------------------------------------------------------------------
		output$driftPlotFixationTime <- renderPlot(
			{
				Plot <- PlotFixationTime(drift_data(), expected_time_flag = input$expected_time_flag)
			}
		)
#------------------------------------------------------------------------------

	}
)
#==============================================================================
