clearvars;

lambda = 380:10:780;

fileName = 'trans_mean_camera.mat';
load(fileName);

% Slicer
vt = reshape(trans_array_m, size(trans_array_m, 1), sizey, sizex);

for wl = 1:size(lambda, 2)
    vvname = sprintf('%s(wl,:,:)','vt');
    vv = squeeze(eval(vvname));
    saveName = insertBefore(fileName, '.mat', ['_' int2str(wl)]);
    save(saveName, 'vv', 'sizey', 'sizex', '-v7.3');
end

% Stacker
clearvars 'sizex' 'sizey';
clearvars 'trans_array_m';

% Load sizex and size y for size of array
loadName = insertBefore(fileName, '.mat', '_1');
load(loadName, 'sizex', 'sizey');
stacked = zeros(size(lambda, 2), sizey, sizex);

for wl = 1:size(lambda, 2)
    loadName = insertBefore(fileName, '.mat', ['_' int2str(wl)]);
    load(loadName);
    stacked(wl, :, :) = vv;
end

trans_array_m = reshape(stacked, 41, sizey * sizex);
save('trans_mean_camera2.mat', 'trans_array_m', 'sizey', 'sizex', '-v7.3');

% Verification: OK
clearvars;
data1 = load('trans_mean_camera.mat');
data2 = load('trans_mean_camera2.mat');
isequaln(data1.trans_array_m, data2.trans_array_m)