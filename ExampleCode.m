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

% Create ReadSpectralDataBase Object
dt = ReadSpectralDatabase;

% Define Standard illuminant
dt.std_d65();

% Set path to BiomaxOrgan 10 data
biomax_path = pwd;
dt.set_biomax_path(biomax_path);

% Set sample name
sample_name = 'Bladder_red';
dt.set_sample_name(sample_name);

% Load the spectral data
dt.load_data;

% Compute CIE coordinates
dt.transmittance2LAB('y') % 'y' to trim the transmittance to 1, can be > 1 in the measurements due to uncertainties

% Compute sRGB values 
disp('Compute sRGB');
dt.XYZ2sRGB;

% Display truth image
disp('Tiff image');
im = dt.img_tiff;
figure;
image(im);
axis image;

% Save the output
% CIE coordinates
XYZ_array = dt.XYZ;
save([pwd '\output\CIE_Coord\XYZ_array'],'XYZ_array');
CovXYZ_array = dt.CovXYZ;
save([pwd '\output\CIE_Coord\CovXYZ_array'],'CovXYZ_array');
LAB_array = dt.LAB;
save([pwd '\output\CIE_Coord\LAB_array'],'LAB_array');
CovLAB_array = dt.CovLAB;
save([pwd '\output\CIE_Coord\CovLAB_array'],'CovLAB_array');

% sRGB and Tiff
rgb = dt.rgb;
save([pwd '\output\RGB\rgb'],'rgb');
truth = im;
save([pwd '\output\RGB\truth'],'truth');