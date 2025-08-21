parameter_data <- tibble(
  Parameter = c(
    'Pairwise Interactions', 'Pairwise Interactions',
    'Pairwise Interactions', 'Pairwise Interactions',
    'HOI Interactions', 'HOI Interactions',
    'HOI Interactions', 'HOI Interactions',
    'Self-regulation'
  ),
  Symbol = c(
    'Mean strength (A_mean)', 'Mean strength (A_mean)',
    'Standard deviation (A_sd)', 'Standard deviation (A_sd)',
    'Mean strength (B_mean)', 'Mean strength (B_mean)',
    'Standard deviation (B_sd)', 'Standard deviation (B_sd)',
    'Intensity (B_diag)'
  ),
  Scenario = c(
    'Grid 1', 'Grid 2',
    'Grid 1', 'Grid 2',
    'Grid 1', 'Grid 2',
    'Grid 1', 'Grid 2',
    'Both Grids'
  ),
  `Min Value` = c(0.05, 0.05, 0.05, 0.1, 0.05, 0.05, 0.05, 0.1, 0),
  `Max Value` = c(1.0, 0.5, 0.1, 0.3, 1.0, 0.5, 0.1, 0.3, 1.5),
  Steps = c(10, 10, 5, 5, 10, 10, 5, 5, 4)
)

param_order <- c('Pairwise Interactions', 'HOI Interactions', 'Self-regulation')
symbol_order <- c('Mean strength (A_mean)', 'Standard deviation (A_sd)',
                  'Mean strength (B_mean)', 'Standard deviation (B_sd)',
                  'Intensity (B_diag)')

parameter_data_new <- parameter_data %>%
  mutate(Parameter = factor(Parameter, levels = param_order),
         Symbol = factor(Symbol, levels = symbol_order))

parameter_table <- parameter_data_new %>%
  gt(groupname_col = 'Parameter', rowname_col = 'Symbol') %>%
  tab_header(title = 'Simulation Parameter Ranges',
             subtitle = 'Parameters were sampled across two distinct grids') %>%
  cols_label(Scenario = 'Scenario',
             `Min Value` = 'Min Value',
             `Max Value` = 'Max Value',
             Steps = 'Steps') %>%
  tab_options(row_group.font.weight = 'bold',
              row_group.background.color = "#f2f2f2")

gtsave(parameter_table, 'parameter_table.png', vwidth = 1000, vheight = 800)