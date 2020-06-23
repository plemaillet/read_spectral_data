# Hyperspectral Microscope (HIMS) tissue transmittance measurements reader

This code was developed to read the hyperspectral transmittance data measured by DIDSR's Hyperspectral Imaging Microscope (HIMS) for a selection of 8 BiomaxOrgan10 tissue microarray slides (US Biomax, 15883 Crabbs Branch Way, MD 20855, USA):

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

This data set will be available soon on the XXX website. The example code outputs the CIE1931 XYZ and CIE1976 L\*a\*b\* coordinates and their covariance matrices. It also outputs the sRGB coordinates of the image.

## Disclaimer

This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.

## Example
"ExampleCode.m" provides an example on how to read and display the input data and the computed outputs for the Uterine Cervix tissue sample. It creates an instance of the ReadSpectralData class that accepts as inputs the hyperspectral data measured by the HIMS (from 380 nm to 780 nm in steps on 10 nm):
- Mean value of the transmittance for each wavelength
- Standard deviation of the transmittance for each wavelength

The outputs for each pixel in the image are:
- The CIE 1931 tri-stimulus coordinates, CIEXYZ
- The covariance matrix on the CIEXYZ coordinates
- The CIEL\*a\*b\* 1976 coordinates
- The covariance matrix on the CIEL\*a\*b\* coordinates
- The sRGB coordinates

The input data are provided for each of the 8 selected tissues in the input/\<Tissue Name\>/Transmittance subfolder. For the example presented here: [input/UterineCervix_red/Transmittance](https://github.com/plemaillet/read_spectral_data/tree/master/input/UterineCervix_red/Transmittance)

The CIEXYZ coordinates of the reference white and the color coordinates (CIEXYZ, CIEYxy and CIEL\*a\*b\*) resulting from the transmittance measurements are displayed for a selection of pixels locations corresponding to superficial cells (X = 266; Y = 52), basal cells (X = 337; Y = 266) and stroma (X = 416; Y = 319). A tiff image is obtained by reshaping the sRGB coordinates:

<p align="center">
  <img width="746" height="512" src="Uterine_Cervix_red_Tagged.png">
</p>

For this selection of pixels positions, the spectral mean values of the transmittance are plotted. A line plot with error bars of the spectral transmittance is also presented (coverage factor *k = 2*, *i.e.* error bars = 2 standard deviations):

<p float="middle">
	<img width="420" height="315" src="plot_T.png" />
	<img width="420" height="315" src="errorbar_T.png" /> 
</p>

The color of these plots can follow the color order from Matlab or the sRGB color computed from the data. The CIEXYZ, CIEL\*a\*b\* and sRGB data can be visualized using a 3D scatter plot with the color of the data point following the sRGB coordinates:

 <p float="middle">
	<img width="420" height="315" src="scatter3_XYZ.png" />
	<img width="420" height="315" src="scatter3_LAB.png" /> 
</p>
 <p align="center">
	<img width="420" height="315" src="scatter3_RGB.png" /> 
</p>

The output CIE coordinates (CIEXYZ, CIEL\*a\*b\* and the respective covariance results) are saved in the [output/UterineCervix_red/CIE_Coord](https://github.com/plemaillet/read_spectral_data/tree/master/output/UterineCervix_red/CIE_Coord) subfolder. The sRGB color coordinates and tiff image are saved in [output/UterineCervix_red/RGB](https://github.com/plemaillet/read_spectral_data/tree/master/output/UterineCervix_red/RGB) subfolder.
