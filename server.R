#==============================================================================
#    server.R : Genetic Drift Simulator Server
#    Copyright (C) 2023  Bruno Toupance <bruno.toupance@mnhn.fr>
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

source("genetic_drift.R")


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
                    simulate_data(
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
        observeEvent(
            input$mul2, 
            {
                updateNumericInput(session, "nb_ind", 
                                   value = input$nb_ind * 2)
                
                updateNumericInput(session, "nb_gen", 
                                   value = input$nb_gen * 2)
            }
        )

#------------------------------------------------------------------------------
        observeEvent(
            input$div2, 
            {
                updateNumericInput(session, "nb_ind", 
                                   value = input$nb_ind %/% 2)
                
                updateNumericInput(session, "nb_gen", 
                                   value = input$nb_gen %/% 2)
            }
        )


#------------------------------------------------------------------------------
        output$export_count <- downloadHandler(
            
                filename = function() {
                    alpha_num <- c(0:9, LETTERS[1:6])
                    rnd_num <- paste(sample(alpha_num, size = 8, 
                                            replace = TRUE), 
                                     collapse = "")
                    file_path <- sprintf("drift_data_%s_%s.txt", Sys.Date(), 
                                         rnd_num)
                    return(file_path)
                },
                
                content = function(file_path) {
                    write.table(
                        get_count_df(drift_data()), 
                        file = file_path, sep = "\t", 
                        quote = FALSE, row.names = FALSE
                    )
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
        output$drift_plot_freq <- renderPlot(
            {
                my_plot <- plot_freq(drift_data(), fix_flag = input$fix_flag)
            }
        )
#------------------------------------------------------------------------------
        output$drift_plot_freq_density <- renderPlot(
            {
                my_plot <- plot_freq_density(
                    drift_data(), 
                    count_flag = input$count_flag)
            }
        )
#------------------------------------------------------------------------------
        output$drift_plot_mean_P <- renderPlot(
            {
                my_plot <- plot_mean_P(
                    drift_data(), 
                    expected_mean_flag = input$expected_mean_flag)
            }
        )
#------------------------------------------------------------------------------
        output$drift_plot_variance_P <- renderPlot(
            {
                my_plot <- plot_variance_P(
                    drift_data(), 
                    expected_var_flag = input$expected_var_flag)
            }
        )
#------------------------------------------------------------------------------
        output$drift_plot_fixation_prob <- renderPlot(
            {
                my_plot <- plot_fixation_prob(
                    drift_data(), 
                    expected_prob_flag = input$expected_prob_flag)
            }
        )
#------------------------------------------------------------------------------
        output$drift_plot_fixation_time <- renderPlot(
            {
                my_plot <- plot_fixation_time(
                    drift_data(), 
                    expected_time_flag = input$expected_time_flag)
            }
        )
#------------------------------------------------------------------------------

    }
)
#==============================================================================
