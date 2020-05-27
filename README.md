# Hyperspectral Microscope (HIMS) tissue transmittance measurements reader

This code was developped to read the transmittance hyperspectral data measured by DIDSR Hyperspectral Imaging MicroScope (HIMS) of a selection of 8 BiomaxOrgan10 tissue microarray slides (US Biomax, 15883 Crabbs Branch Way, MD 20855, USA):

| Tissue Name | Position |
| --- | :---: |
| Bladder | M13 |
| Brain | H10 |
| Breast | A1 |
| Colon | H6 |
| Kidney | H7 |
| Liver | H9 |
| Lung | J7 |
| Uterine Cervix | B10 |

It outputs the CIE1931 XYZ and CIE1976 L*a*b* coordinates and their covariance matrices. It also outputs the sRGB coordinates and the tissue color truth as a tiff image.

## Disclaimer

This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.

## Example simulation

Example Bladder (image reduced in size to 557 x 512 pixels

Simulation input files, and auxiliary files, have been added in the [example_simulations folder](https://github.com/DIDSR/VICTRE_MCGPU/tree/master/example_simulations) to allow the replication of two of the simulations of the VICTRE project.
A breast phantom with scattered glandularity and a heterogreously dense phantom created with C. Graff model are provided (these phantoms were not part of the original VICTRE pivotal study, and have a single large mass embedded inside). The material files (generated from PENELOPE 2006 material files) and energy spectra used in the simulations are included.

<img src="output/Bladder_red/RGB/truth.tif" width="557px"/>

### (e) Some list
List:

  * Item 1.
  * Item 2.
