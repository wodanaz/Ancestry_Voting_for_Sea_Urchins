First, we need to create a conda environmet in your HPS system. This means, Anaconda should be supported.


``` bash
module load Anaconda3

conda env create -f Seurat4.yml --prefix /data/wraycompute/alejo/conda/Seurat4 

```


note you need a Seurat4.yml file, to create simply run:

```bash 
nano Seurat4.yml
channels:
  - conda-forge
  - bioconda
  - r
dependencies:
  - r-base
  - r-seurat
  - python=3.6
  - r-monocle3
  - r-fnn
  - r-tidyverse
  - r-lattice
  - r-viridis
  - r-reshape
  - r-doparallel
  - r-foreach
  - r-ape
  - agat
```

Then save it and exit
