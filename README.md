# CellRomeR

This is an R package for phenotypic clustering of cells based on their morphological and motility features from microscopy images. This takes live cell microscopy imaging data after segmentation and tracking of the cells as input. The current implementation is compatible with morphological and motility data exported from the TrackMate tool as xml files.


**Installation instructions**
``` R 
install.packages("devtools")
devtools::install_github("elolab/cellromer")

```

Or from source (in your shell)
```
git clone THISPAGE
R -e 'install.packages("devtools"); devtools::install("cellromer")'
```


