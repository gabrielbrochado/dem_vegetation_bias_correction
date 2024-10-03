# DEM Vegetation Bias Correction and Comparison Scripts

This is the repository of the python scripts used in the article "New Method to Correct Vegetation Bias in Copernicus Digital Elevation Model to Improve Flow Path Delineation" by Gabriel Thomé Brochado and Camilo Daleles Rennó.

The input and output data of these scripts were made available as suplementary materials (zip file), and can be downloaded at: https://drive.google.com/file/d/1ryoEhMyneMQtiplpiE6Y1CQdBA2a8Uqp/view?usp=sharing

To reproduce the results, this file should be unziped in a base folder and the scripts executed in the given order. However, the script used to generate the reference flow paths requires as input the vectorial drainage networks of the study areas, which were not made available with the suplmentary materias, because its distribution requires the authorization of the Brazilian Army Geographic Service. Alternatevely, the reference flow paths present in the suplementary materials can be used as inputs for the next step.

The first script is the DEM correction method it self, while the others constitute the comparison method used to evalute the performance of the original Copernicus DEM, the benchmark (FABDEM) and the new corrected DEM, regarding stream flow path delineation.

## Scripts Description

### 1-dem correction
Corrects the vegetation bias on Copernicus DEM.

**Input datasets**: Copernicus DEM, forest height, forest loss per year.

**Outputs**: Corrected Copernicus DEM and folder with various auxiliary files.

### 2-generate ref flowpaths
Generates the set of reference flow paths based on the vectorial drainage networks of the study areas. These flow paths are used on the DEM comparison method applied in the study.

**Input**: drainage network shapefiles of all study areas, and the adjusted forest height dataset.

**Outputs**: reference flow paths for each area and radius.

### 3-generate dem displacement area polygons
Generates the displacement area polygons for each reference flow path and each DEM, for all areas and radius.

**Inputs**: reference flow paths for each area and radius.

**Outputs**: displacement area polygons for each area and radius, and a dataframe compiling the results.

### 4-analyze displacement areas
Compares the displacement areas of the DEMs, for each area and radius, using the wilcoxon signed-rank test.

**Input**: dataframe with the displacement areas of the DEMs for each area and radius.

**Outputs**: dataframe with the wilcoxon test p-values resulting from the comparison of each pair of DEMs for each area and radius; filtered dataframe with the samples used on the tests.

### 5-generate sampled displacement area polygons
Filters the set of displacement area polygons of each area and radius by preserving only the ones used on the tests.

**Input**: filtered dataframe with the displacement areas used on the tests.

**Output**: displacement area polygons for each area and radius, that were used in the tests.





