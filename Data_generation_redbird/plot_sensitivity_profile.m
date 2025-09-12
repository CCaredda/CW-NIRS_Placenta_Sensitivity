
clear
close all
clc

addpath('~/Soft/redbird/matlab');
addpath('~/Soft/iso2mesh');
addpath('../functions');

%Display results
display = 0;

%Wavelength in nm
Lambdas = 780;

%Volume ize in mm
xdim_mm = 200;
ydim_mm = 200;
zdim_mm = 200;
% max_vol_mesh = [0.1; 0.1; 1; 1000];
max_vol_mesh = [0.5; 0.5; 1; 10000];

%Thickness layer
thickness_layers_mm = [2 4 10];

%Source detector separation in mm
detectors_SD_mm = [30 40 50];

clear cfg;

%Create 4 layers volume
cfg = create_meshed_volume_4layers(xdim_mm, ydim_mm, zdim_mm, thickness_layers_mm, max_vol_mesh, detectors_SD_mm, display);

C_HbT_muscle = 25*1e-6;
C_HbT_placenta = 25*1e-6;
SatO2_muscle = 0.6;
SatO2_placenta = 0.8;
% f_mel = 0.0255;
f_mel = 0.305;

% Calculate optical properties
optical_prop = process_optical_properties_skin_Fat_muscle_placenta(Lambdas,f_mel,SatO2_muscle, SatO2_placenta,C_HbT_muscle,C_HbT_placenta);




%Calculate sensisitivity profile
[sensitivity_profile, Diffuse_reflectance] = get_sensitivity_profiles(cfg, optical_prop);

%Calculate sensiticity indexes
% S_index = get_sensitivity_index(cfg, sensitivity_profile, thickness_layers_mm);


%Init output
Sensitivity_proba = zeros(xdim_mm,ydim_mm,zdim_mm,size(sensitivity_profile,2));

%interpolate volume
for d=1:size(sensitivity_profile,2)
    [xi, yi, zi] = meshgrid(0.5:xdim_mm-0.5, 0.5:ydim_mm-0.5, 0.5:zdim_mm-0.5);
    Sensitivity_proba(:,:,:,d) = griddata(cfg.node(:,1), cfg.node(:,2), cfg.node(:,3), sensitivity_profile(:,d), xi, yi, zi);

end
close all
figure(1)
subplot(131)
imagesc(log10(Sensitivity_proba(:,:,1,1)))
subplot(132)
imagesc(log10(squeeze(Sensitivity_proba(:,100,:,1))))
subplot(133)
imagesc(log10(squeeze(Sensitivity_proba(100,:,:,1))))

srcpos = cfg.srcpos;
detpos = cfg.detpos;

save(strcat("Sensitivity_proba_fmel_",num2str(f_mel),".mat"),'Sensitivity_proba', 'srcpos', 'detpos', 'thickness_layers_mm');