%%%%%%%
% Example code to load the data and output the CIE coordiantes and sRGB results
%%%

% input = 
%   - transmittance:                                    41 x (sizey*sizex) array
%   - standard deviation on transmittance:              41 x (sizey*sizex) array
% output = 
%   - CIEXYZ: CIE 1931 tri-stimulus coordinates         (sizey*sizex) x 3 array
%   - CovXYZ: Covariance matrix on XYZ coordinates      3 x 3 x (sizey*sizex) array
%   - CIELAB: CIELab 1976 coordinates                   (sizey*sizex) x 3 array
%   - CovLAB: Covraince matrix on Lab coordinates       3 x 3 x (sizey*sizex) array
%   - rgb: sRGB                                         (sizey*sizex) x 3 array
%   - im                                                sizey x sizex x 3 array

clearvars;
close all;

% Path to BiomaxOrgan 10 data
biomax_path = pwd;

% Name of the sample
sample_name = 'Bladder_red';

% Create ReadSpectralDataBase Object
dt = ReadSpectralDatabase(biomax_path, sample_name);

% Define Standard illuminant
d65 = LightSource([pwd '\input\DataIlluminants\spec_cied65']);
dt.set_ls(d65.ls);

% Compute CIEXYZ coordinates
dt.transmittance2XYZ('y') % 'y' to trim the transmittance to 1, can be > 1 in the measurements due to uncertainties

% Compute CIELAB coordinates
dt.transmittance2LAB('y') % 'y' to trim the transmittance to 1

% Compute sRGB values 
dt.transmittance2sRGB ('y') % 'y' to trim the transmittance to 1

% Reshape sRGB to tiff and display truth image
disp('Tiff image');
im = dt.img_tiff;
figure;
image(im);
axis image;

% Save the outputs
% CIE coordinates
XYZ_array = dt.XYZ;
save([pwd '\output\Bladder_red\CIE_Coord\XYZ_array'],'XYZ_array');
CovXYZ_array = dt.CovXYZ;
save([pwd '\output\Bladder_red\CIE_Coord\CovXYZ_array'],'CovXYZ_array');
LAB_array = dt.LAB;
save([pwd '\output\Bladder_red\CIE_Coord\LAB_array'],'LAB_array');
CovLAB_array = dt.CovLAB;
save([pwd '\output\Bladder_red\CIE_Coord\CovLAB_array'],'CovLAB_array');

% sRGB and Tiff
rgb = dt.rgb;
save([pwd '\output\Bladder_red\RGB\rgb'],'rgb');
imwrite(im ,[pwd '\output\Bladder_red\RGB\truth.tif']);