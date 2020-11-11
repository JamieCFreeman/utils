
# 2020-11-11 JCF

# Would like to start saving color palette library across projects, b/c choosing
#   colors is very time consuming.

# Will use this as a master file to save color palettes I use in publication figures.
#   Can source this file in other scripts and use palettes by name.

# Potential improvement: Maybe a function that would print all palette names or 
#   plot all palettes


# Function for visualizing palettes, for a given palette, plots image of all colors.
see_palette <- function(pal) {
  
  image(1:length(pal), 1, as.matrix(1:length(pal)), 
        col=pal,
        xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n" , main=substitute(pal))
  
  
  }


# Color palette generator I like:
# https://coolors.co/801515-d3d0cb-e2c044-587b7f-1d2f6f

# Color palette used for SV fig
dark.primary <-  c("#801515", "#D3D0CB", "#E2C044", "#587B7F", "#1D2F6F")


# Color palette used for BSA strains
# KS17-R. IsoCS, 
BSA.strains <-  c("#FFA915", "#2A7FFF", "#55286F", "#59C9A5", "#28231C")




see_palette(dark.primary)


see_palette(BSA.strains)
