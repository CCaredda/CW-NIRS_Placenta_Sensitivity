clear
close all
clc


run_in_cluster = 1;

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

%Wavelength
Lambdas = 780;

%Arrays thickness
skin_thickness = 1;
adipose_thickness = 4;
muscle_thickness = 10;

%Saturation array
SatO2_muscle = 0.6;
SatO2_placenta = 0.8;


%HbT array
C_HbT_muscle = 50*1e-6;
C_HbT_placenta = 50*1e-6;

%Volume fraction of melanosome according to the color tones
% Modeling and Verification of Melanin Concentration on Human Skin Type
f_mel = 0.0255;


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

% Add path
addpath('functions');
if run_in_cluster == 1
    addpath('/pbs/home/c/ccaredda/private/mcx/utils');
    addpath('/pbs/home/c/ccaredda/private/mcxlab');
else
    addpath('/home/caredda/Soft/mcx/utils');
    addpath('/home/caredda/Soft/mcxlab');
    % addpath('/home/caredda/Soft/iso2mesh');
    % addpath('/home/admin/Software/mcx-linux-x86_64-v2023/mcx/utils');
    % addpath('/home/admin/Software/mcxlab-linux-x86_64-v2023/mcxlab');  
end


% Source from 0
cfg.issrcfrom0=1;

%Volume square size
volume_square_size = 200;

% Source type (values in accordance with VOlume size, se below)
cfg.srctype='pencil';
src_pos = [floor(volume_square_size/2 -1) floor((volume_square_size - 80/model_resolution_in_mm)/2 -1) 0];
cfg.srcdir=[0 0 1];

%Define multiple detectors  (values in accordance with VOlume size, se below)
det_pos = zeros(8,4);
for i=1:8
    det_pos(i,:) = [src_pos(1) src_pos(2)+i*floor(10/model_resolution_in_mm) 1 0.5]; % radius: 0.5 mm
end
cfg.detpos = det_pos;

cfg.savedetflag = 'dp'; %Save detector id and partial path length

% GPU processing
cfg.gpuid=1;


%Save constant
if ~exist('output_lookup_table_skin', 'dir')
    mkdir('output_lookup_table_skin')
end

% Use the dir function to list files based on the pattern
%fileList = dir(strcat('output_lookup_table/out_St_muscle_',num2str(SatO2_muscle),'_St_placenta_',num2str(SatO2_placenta),'*_HbT_muscle_umol_',num2str(C_HbT_muscle*1e6),'_HbT_placenta_umol_',num2str(C_HbT_placenta*1e6),'.mat'));

% Count the number of files
%skip_items = numel(fileList);




                    
% Define optical properties
optical_prop = process_optical_properties_skin_Fat_muscle_placenta(Lambdas,f_mel,SatO2_muscle, SatO2_placenta,C_HbT_muscle,C_HbT_placenta);

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
N_skin = floor(skin_thickness/model_resolution_in_mm);
N_adipose = floor(adipose_thickness/model_resolution_in_mm);
N_muscle = floor(muscle_thickness/model_resolution_in_mm);

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

nb_noise = 1;
array_Diffuse_reflectance = zeros(nb_noise,8);
array_Sensitivity_index = zeros(nb_noise,8,4);


for id_noise = 1:nb_noise

	%Process simulation with light source at the source
	%position
	cfg.srcpos=src_pos;

	% Random seed to obtain different results when running multiple simulations for the same input parameters
	cfg.seed = randi([0,99999],1);

	% calculate the fluence and partial path lengths
	[fluence,output_det]=mcxlab(cfg);
	% ppath = output_det.ppath;
	% detid = output_det.detid;
	  
	% Detector pos
	for i=1:length(det_pos)
	    %Diffuse reflectance
	    array_Diffuse_reflectance(id_noise,i) = fluence.dref(det_pos(i,1)+1,det_pos(i,2)+1,1);
	end

	%Add numbers of detected photons


	% Compute simulation at detector location
	for d = 1:size(det_pos,1)
	    % Random seed to obtain different results when running multiple simulations for the same input parameters
	    cfg.seed = randi([0,99999],1);
	    
	    %Set light source position
	    cfg.srcpos= [src_pos(1) det_pos(d,2) src_pos(3)];

	    % calculate the fluence
	    [flux]=mcxlab(cfg);

	    %Calculate sensitivity profile
	    Sensitivity_profile = fluence.data.*flux.data;

	    %Normalize by the sum to get the density
	    %probability
	    Sensitivity_profile = Sensitivity_profile/sum(Sensitivity_profile,"all");
        Sensitivity_profile(isnan(Sensitivity_profile)) = 0;

        %Skin
	    array_Sensitivity_index(id_noise,d,1) = sum(Sensitivity_profile(:,:,2:2+N_skin-1),"all");	
        %Adipose tissue
    	array_Sensitivity_index(id_noise,d,2) = sum(Sensitivity_profile(:,:,2+N_skin:2+N_adipose+N_skin-1),"all");
        %Muscle
        array_Sensitivity_index(id_noise,d,3) = sum(Sensitivity_profile(:,:,2+N_skin+N_adipose:2+N_skin+N_adipose+N_muscle-1),"all");
        %Placenta
        array_Sensitivity_index(id_noise,d,4) = sum(Sensitivity_profile(:,:,2+N_skin+N_adipose+N_muscle:end),"all");
	end
end

save('output_noise.mat','array_Diffuse_reflectance','array_Sensitivity_index');                

                            









