clear
close all
clc

cluster = 0;

% Add path
if cluster ==1
    addpath('~/private/redbird/matlab');
    addpath('~/private/iso2mesh');
else
    addpath('~/Soft/redbird/matlab');
    addpath('~/Soft/iso2mesh');
end

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

%Source detector separation in mm
detectors_SD_mm = [30 40 50];


%Arrays thickness
thickness_layers_mm = [4 2 3];

%optical prop
optical_prop.mua_skin = 0.1*0.05; % mm-1
optical_prop.mu_s_skin =  0.1*9.3; % mm-1
optical_prop.g_skin= 0.9;
optical_prop.n_skin =1.37;

optical_prop.mua_adipose = 0.1*0.98; % mm-1
optical_prop.mu_s_adipose = 0.1*12; % mm-1
optical_prop.g_adipose= 0.9;
optical_prop.n_adipose=1.37;

optical_prop.mua_muscle = 0.1*0.11; % mm-1
optical_prop.mu_s_muscle= 0.1*8; % mm-1
optical_prop.g_muscle= 0.9;
optical_prop.n_muscle=1.37;
    
optical_prop.mua_placenta = 0.1*0.1; % mm-1
optical_prop.mu_s_placenta=  0.1*10; % mm-1
optical_prop.g_placenta =0.9;
optical_prop.n_placenta=1.37;
    



%Create 4 layers volume
cfg = create_meshed_volume_4layers(xdim_mm, ydim_mm, zdim_mm, thickness_layers_mm, max_vol_mesh, detectors_SD_mm, display);

%Get Sensitivity indexes for the 4 layers
fprintf(1,strcat('Calculating sensitivity index\n'));

%Calculate sensisitivity profile
[sensitivity_profile, Diffuse_reflectance] = get_sensitivity_profiles(cfg, optical_prop);

%Calculate sensiticity indexes
Sensitivity_indexes = get_sensitivity_index(cfg, sensitivity_profile, thickness_layers_mm);


save(strcat('Phantom_Data_',num2str(Lambdas),'.mat'),'Diffuse_reflectance','Sensitivity_indexes');


