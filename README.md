# CellRomeR

This is an R package for phenotypic clustering of cells based on their morphological and motility features from microscopy images. This takes live cell microscopy imaging data after segmentation and tracking of the cells as input. The current implementation is compatible with morphological and motility data exported from the TrackMate tool as xml files.

CellRomeR is an R package for the analysis and visualisation of cell tracking data from TrackMate [1]. It uses ILoReg [2] for clustering and includes several different plotting functions. Input data is an XML file from TrackMate.

#### Installation
``` R 
install.packages("devtools")
devtools::install_github("elolab/cellromer")

```

Or from source (in your shell)
```
git clone THISPAGE
R -e 'install.packages("devtools"); devtools::install("cellromer")'
```


#### Example analysis

Load the library
``` R
library(CellRomeR)
```

##### Data Import

XML is the native data format used by TrackMate to store data and the models. In TrackMate data, long tracks can often branch one or more times due to cell divisions or tracking issues. Plotting these branching tracks can pose challenges. To address this, it is advisable to only export non-branching tracks from TrackMate. We provide the `IRmigr.TMxml()` function to read single TrackMate data into `MigrDat` migration data object. It reads and stores spots, tracks, and edges data in separate slots in a `MigrDat` object.

An example XML output file form TrackMate is available here: LINK
The example data is based on the cancer cell migration dataset available in [Zenodo](https://zenodo.org/records/5206107). The original trackMate `P31-crop.xml` dataset was imported into TrackMate (load in TrackMate file), re-tracked, and saved with as non-branching version (`P31-crop-nb.xml`).


``` R
migrdata <- import_XML("P31-crop_nb.xml")
migrdata <- Apply_TMate_Filter(migrdata)
```

##### `MigrDat` data object

The data is stored in a MigrDat object which has the following slots:

- Spots : Containing information related to each spot
- Tracks : Information regarding tracks and movement of the cells
- Edges : Information regarding edges ## Iivari please check!
- Roi_points : Defining *outlines of the spots/cells* (regions of interest)
- Metadata : Containing any extra information that might be useful
- IDs : The provided IDs from the TrackMate output
- Dimensional reductions: Containing calculated distance coordinates of spots, tracks, and edges
    
##### Clustering

The function `clustering()` can be used to run dimensionality reduction and clustering. The clustering is based on ILoReg [2], an ensembl clustering method developed for clustering of single-cell RNA-sequencing data. The argument `uniq` is used to store the clustering result with a unique identifier for later identification. NOTE! Clustering may take several minutes depending on the size of the data and computational resources.

``` R
migrdata <- clustering(MigrObj=migrdata, dat.slot="raw", type="S", vars=features, kILoReg=5)
```

The clustering above was done for Spot (S) type only, but by adjusting the `type` parameter we can also do this for Tracks (T) and Edges (E).

##### Visualisation

cellRomeR includes several plotting options. The `plot` function plots the trajectories of the cells. The `plot_violin` can be used on any of the data features. 

``` R
plot(migrdata)
plot_umap(migrdata)
plot_pie(migrdata)
plot_heatmap(migrdata)
plot_violin(migrdata, feature="CIRCULARITY")
plot_rag(migrdata)
```

#### Citations

1. Ershov, D., Phan, M.-S., Pylvänäinen, J. W., Rigaud, S. U., et al (2022). TrackMate 7: integrating state-of-the-art segmentation algorithms into tracking pipelines. Nature Methods, 19(7), 829–832. doi:10.1038/s41592-022-01507-1

2. Smolander J., et al (2021). ILoReg: a tool for high-resolution cell population identification from single-cell RNA-seq data. Bioinformatics, 37(8):1107-1114. doi: 10.1093/bioinformatics/btaa919.




