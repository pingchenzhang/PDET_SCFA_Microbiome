
draw_venn_diagram <- function(
  a = 'a',
  b = 'b',
  ab = c('a', 'b'),
  a_name = 'A',
  b_name = 'B',
  main = 'Venn',
  fsize = 8, fsize_main=10
) {
  # Define the circle's center position and radius
  a_center <- c(x = 3, y = 5)
  b_center <- c(x = 7, y = 5)
  radius <- 5
  df <- data.frame()
  angles <- seq(0, 2 * pi, length.out = 100)
  ggplot(df, aes(x = x, y = y)) +
    geom_path(
      data = data.frame(
        x = a_center["x"] + radius * cos(angles),
        y = a_center["y"] + radius * sin(angles)
      ),
      aes(x = x, y = y),
      size = 1
    ) +
    geom_path(
      data = data.frame(
        x = b_center["x"] + radius * cos(angles),
        y = b_center["y"] + radius * sin(angles)
      ),
      aes(x = x, y = y),
      size = 1
    ) +
    # text label
    annotate("text", x = 0, y = 5, label = paste(a, collapse = '\n'), size = fsize) +
    annotate("text", x = 5, y = 5, label = paste(ab, collapse = '\n'), size = fsize) +
    annotate("text", x = 10, y = 5, label = paste(b, collapse = '\n'), size = fsize) +
    annotate("text", x = 2, y = -.5, label = a_name, size = fsize + 1) +
    annotate("text", x = 8, y = -.5, label = b_name, size = fsize + 1) +
    annotate("text", x = 5, y = 11, label = main, size = fsize_main) +
    # Graphic settings
    #xlim(0, 10) +
    #ylim(0, 10) +
    theme_void() +
    theme(
      # plot.margin = margin(20, 20, 20, 20),
      panel.background = element_rect(fill = "white", color = NA)
    )
}
# draw_venn_diagram()
