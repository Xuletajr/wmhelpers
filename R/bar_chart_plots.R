# =============================================================================
#            Function to automate plots and avoid duplicate code 
# =============================================================================

# No need to edit. Only if necessary an update

# var1 = count of interest
# var2 =  species name.
# var3 = IUCN status
# Var4 = Species of greatest conservation need. Specific for some project

bar_chart <- function(dat, var1, var2, var3, var4 = NULL, xlab, subtitle,
                      breaks_x = NULL, limits_x = NULL) {
  
  # Create a new label column in the data frame based on the conditions for var3 and var4
  dat <- dat %>%
    dplyr::mutate(
      # Check conditions for var3 and var4, combine them accordingly
      label = case_when(
        {{ var3 }} == "LC" ~ ifelse({{ var4 }} == "Yes", "*", ""),  # If var3 is "LC", just show '*' if var4 is "Yes"
        {{ var4 }} == "Yes" ~ paste0({{ var3 }}, "*"),  # If var4 is "Yes", concatenate var3 and "*"
        TRUE ~ as.character({{ var3 }})  # Otherwise, just use var3
      )
    )
  
  # Create the base plot
  p <- dat %>%
    ggplot2::ggplot(aes(
      x = {{ var2 }},
      y = reorder({{ var1 }}, {{ var2 }}))) +
    
    geom_bar(
      stat = "identity",
      fill = "#00D973") +
    
    # Plot the labels using the newly created 'label' column
    geom_text(aes(label = label), hjust = -0.2, size = 4) +
    
    labs(
      x = xlab,
      y = NULL,
      subtitle = subtitle) +
    theme_bw() +
    theme(
      axis.title = element_text(
        size = 18,
        color = config$color2),
      axis.text = element_text(
        color = config$color2,
        size = 12),
      plot.subtitle = element_text(
        size = 16,
        color = config$color2,
        hjust = 1.0,
        face = "bold"),
      legend.position = "None",
      plot.margin = margin(1, 0.5, 1, 0.5, "cm"))
  
  # Sometimes the x-axis does not come out with good formatting
  # so we can change with the number of breaks
  if (!is.null(breaks_x)) {
    p <- p + scale_x_continuous(
      labels = label_number(accuracy = 1),
      n.breaks = breaks_x,
      limits = limits_x)
  } else {
    p <- p + scale_x_continuous(
      labels = label_number(accuracy = 1),
      expand = expansion(c(0.01, 0.1)))
  }
  
  return(p)
}
