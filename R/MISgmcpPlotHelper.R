#---------------- -
# Function to generate color codes
#---------------- -
generate_colors <- function(n, palette_name = "Pastel2") {
  num_colors <- max(3, n)
  color_palette <- RColorBrewer::brewer.pal(n = num_colors, name = palette_name)
  return(color_palette[1:n])
}
#---------------- -
