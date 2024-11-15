my_plot <- function(data) {
  
  if("L1" %in% colnames(data)){
    plotdata <- rbind(data, data.frame("ID" = unique(data$ID), L0 = unique(data$L0),
                                       "Time" = 0, "Delta" = "start", L1 = 0, A = unique(data$A)))
  }                                  
  else{
    plotdata <- rbind(data, data.frame("ID" = unique(data$ID), L0 = unique(data$L0),
                                       "Time" = 0, "Delta" = "start"))
  }
  
  diff_events <- plotdata$Delta %>% unique() %>% length()
  cols <- c("darkorange", "blue3", "lightgreen", "darkred", "pink")
  shapess <- c(8, 20, 20 , 8, 20)
  
  ggplot(plotdata) +
    geom_line(aes(x = Time, y = ID, group = ID), color = "grey60", size = 0.7) +  # Light gray lines for each ID
    geom_point(aes(x = Time, y = ID, shape = factor(Delta), color = factor(Delta)), size = 3) + # Use shapes and colors for Delta
    theme_minimal(base_size = 15) +  # Increase base font size for readability
    scale_shape_manual(values = shapess[1:diff_events]) +  # Customize shapes for different events
    scale_color_manual(values = cols[1:diff_events]) +  # Customize colors for events
    labs(
      title = "Survival Data",
      x = "Time",
      y = "Patient ID",
      shape = "Event Type",
      color = "Event Type"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold title
      axis.title.y = element_text(margin = margin(r = 10)),  # Add space to y-axis title
      axis.title.x = element_text(margin = margin(t = 10)),  # Add space to x-axis title
      legend.position = "top"  # Place legend on top for easy viewing
    )
}