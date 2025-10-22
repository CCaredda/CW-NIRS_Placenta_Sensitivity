clear
close all
clc

run_in_cluster = 0;

% Add path
addpath('../functions');
if run_in_cluster == 1
    addpath('/pbs/home/c/ccaredda/private/mcx/utils');
    addpath('/pbs/home/c/ccaredda/private/mcxlab');
    addpath('~/private/redbird/matlab');
    addpath('~/private/iso2mesh');
else
    addpath('~/Soft/mcx/utils');
    addpath('~/Soft/mcxlab');
    addpath('~/Soft/iso2mesh');
    addpath('~/Soft/redbird/matlab');
    addpath('~/Soft/iso2mesh');
end

%Wavelength
Lambdas = 780;

%Arrays thickness
thickness_skin_in_mm = 3;
thickness_adipose_in_mm = 5;
thickness_muscle_in_mm = 12;

%Saturation array
SatO2_muscle = 0.6;
SatO2_placenta = 0.8;

%HbT array
C_HbT_muscle = 25*1e-6;
C_HbT_placenta = 35*1e-6;


%Volume fraction of melanosome according to the color tones
% Modeling and Verification of Melanin Concentration on Human Skin Type
f_melanosome = 0.155;


%Volume square size
volume_square_size = 200;

%Volume ize in mm
xdim_mm = volume_square_size;
ydim_mm = volume_square_size;
zdim_mm = volume_square_size;
% max_vol_mesh = [0.1; 0.1; 1; 1000];
max_vol_mesh = [0.5; 0.5; 1; 10000];

%Info detectors
detectors_SD_mm = [30, 40, 50];

% Define optical properties
optical_prop = process_optical_properties_skin_Fat_muscle_placenta(Lambdas,f_melanosome,SatO2_muscle, SatO2_placenta,C_HbT_muscle,C_HbT_placenta);


%Detector fiber optics
radius_fiber_optics_mm = 2.3;

out_dir = 'compare_MCX_Redbird';
if ~exist(out_dir, 'dir')
    mkdir(out_dir)
end


%% MCX SIMULATIONS %%

% Repeat the simulation x times
repetitions = 1;
cfg.respin = repetitions;

% Number of photons
nphotons = 1e9;
cfg.nphoton = nphotons; 

% Maximum number of photons that can be detected 
cfg.maxdetphoton = 1e8;

%Size of 1 voxel in mm
model_resolution_in_mm = 1;

% Save diffuse reflectance
cfg.issaveref = 1;
                                           
%Acquisition time
cfg.tstart=0; % Starting time of the simulation (in seconds)
cfg.tend=1e-9; % Ending time of the simulation (in seconds)
cfg.tstep=1e-9; % Time-gate width of the simulation (in seconds)

% Calculate specular reflection if source is outside        
cfg.isspecular = 1; 
cfg.autopilot = 1;
   
% Voxel size in mm
cfg.unitinmm = model_resolution_in_mm;


% Source from 0
cfg.issrcfrom0=1;


%Light source
src_pos = [(xdim_mm/2-1) (ydim_mm - detectors_SD_mm(end))/2-1 0];
cfg.srcdir=[0 0 1];
cfg.srctype = 'pencil';

%Detectors
det_pos = zeros(length(detectors_SD_mm),4);
for i=1:length(detectors_SD_mm)
    det_pos(i,:) = [src_pos(1) src_pos(2)+detectors_SD_mm(i) 1 radius_fiber_optics_mm];
end
cfg.detpos = det_pos;


cfg.savedetflag = 'dp'; %Save detector id and partial path length

% GPU processing
cfg.gpuid=1;


% Set optical properties % [mua,mus,g,n]
% 0: Air
% 1: Skin
% 2: Adipose tissue
% 3: Muscle
% 4: Placenta
cfg.prop=[0 0 1 1; ...
optical_prop.mua_skin optical_prop.mu_s_skin optical_prop.g_skin optical_prop.n_skin ; ...
optical_prop.mua_adipose optical_prop.mu_s_adipose optical_prop.g_adipose optical_prop.n_adipose ; ...
optical_prop.mua_muscle optical_prop.mu_s_muscle optical_prop.g_muscle optical_prop.n_muscle ; ...
optical_prop.mua_placenta optical_prop.mu_s_placenta optical_prop.g_placenta optical_prop.n_placenta];

        

            
%Indexes of layers along z axis
N_skin = floor(thickness_skin_in_mm/model_resolution_in_mm);
N_adipose = floor(thickness_adipose_in_mm/model_resolution_in_mm);
N_muscle = floor(thickness_muscle_in_mm/model_resolution_in_mm);

% Create volume
cfg.vol = 4*ones(volume_square_size,volume_square_size,volume_square_size); %200 mm x 200 mm x 200 mm

% 0: Air
cfg.vol(:,:,1) = 0;
% 1: skin
cfg.vol(:,:,2                   :2+N_skin-1) = 1; % Skin
% 2: Adipose tissue
cfg.vol(:,:,2+N_skin            :2+N_adipose+N_skin-1) = 2; % Adipose
% 3: Muscle
cfg.vol(:,:,2+N_skin+N_adipose  :2+N_skin+N_adipose+N_muscle-1) = 3; % Muscle
% 4: rest is Placenta

%Convert in uint8
cfg.vol = uint8(cfg.vol);

%Slice of the tissue
slice_tissue = squeeze(cfg.vol(floor(volume_square_size/2),floor(volume_square_size/2),:));


%Init output
Diffuse_reflectance = zeros(size(det_pos,1),1);
Sensitivity_indexes = zeros(size(det_pos,1),4);
Sensitivity_profile = zeros(xdim_mm,ydim_mm,zdim_mm,size(det_pos,1));


%Process simulation with light source at the source
%position
cfg.srcpos=src_pos;

% Random seed to obtain different results when running multiple simulations for the same input parameters
cfg.seed = randi([0,99999],1);

% calculate the fluence and partial path lengths
[fluence,output_det]=mcxlab(cfg);
  
% Detector pos
for i=1:length(det_pos)
    %Diffuse reflectance
    Diffuse_reflectance(i) = fluence.dref(det_pos(i,1)+1,det_pos(i,2)+1,1);
end

% Compute simulation at detector location
for d = 1:size(det_pos,1)
    % Random seed to obtain different results when running multiple simulations for the same input parameters
    cfg.seed = randi([0,99999],1);
    
    %Set light source position
    cfg.srcpos= [src_pos(1) det_pos(d,2) src_pos(3)];

    % calculate the fluence
    [flux]=mcxlab(cfg);

    %Calculate sensitivity profile
    S = fluence.data.*flux.data;

    %Normalize by the sum to get the density
    %probability
    S = S/sum(S,"all");
    Sensitivity_profile(:,:,:,d) = S;


    %Get sensitivity indexes
    for m=1:4
        % Find indices of tissue m
        indices = find(slice_tissue == m);
        % Compute sensitivity
        Sensitivity_indexes(d,m) = sum(S(:, :, indices),"all");
    end
end

%Save outputs
output_name = strcat(out_dir,'/MCX.mat');
save(output_name,'Diffuse_reflectance','Sensitivity_indexes','Sensitivity_profile','src_pos','det_pos');


%% PROCESS REDBIRD %%

display = false;




clear cfg
%Thickness layer
thickness_layers_mm = [thickness_skin_in_mm thickness_adipose_in_mm thickness_muscle_in_mm];
    
%Create 4 layers volume
cfg = create_meshed_volume_4layers(xdim_mm, ydim_mm, zdim_mm, thickness_layers_mm, max_vol_mesh, detectors_SD_mm, display);

%Model fiber detector
cfg.radius_fiber_det_mm = radius_fiber_optics_mm;
cfg.reso_detector_mm = 0.1;


%Calculate sensisitivity profile
[sensitivity_profile, Phi_at_det, DR_at_fiber_detector] = get_sensitivity_profiles(cfg, optical_prop);

%Calculate sensiticity indexes
Sensitivity_indexes = get_sensitivity_index(cfg, sensitivity_profile, thickness_layers_mm);

%Save outputs
src_pos = cfg.srcpos;
det_pos = cfg.detpos;
output_name = strcat(out_dir,'/Redbird.mat');
save(output_name,'Phi_at_det','Sensitivity_indexes', 'DR_at_fiber_detector', 'sensitivity_profile','src_pos','det_pos');