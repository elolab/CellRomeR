plot_roi <- function(MigrObj, idx,spots.slot="raw")
  for(roiX in MigrObj@roi_points[[spots.slot]][idx]){
    x_coords <- roiX[c(TRUE,FALSE)]
    y_coords <- roiX[c(FALSE,TRUE)]
    plot(x = x_coords, y = y_coords)
    points(x = x_coords, y = y_coords, pch = 21, bg = "red", col = "black")  
    lines(x = c(x_coords, x_coords[1]), y = c(y_coords, y_coords[1]), col = "blue" )  
  }

