% TO DO calculate bananashape 
% get diffuse reflectance
clear
close all

run_in_cluster = 0;
Nb_measure = 1000;

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
L1_thickness_mm = 2;
L2_thickness_mm = 4;
L3_thickness_mm = 3;


%optical prop
mua_L1 = 0.1*0.98; % mm-1
mua_L2 = 0.1*0.05; % mm-1
mua_L3 = 0.1*0.11; % mm-1
mua_Bulk = 0.1*0.9; % mm-1

mus_L1 = 0.1*12; % mm-1
mus_L2 = 0.1*9.3; % mm-1
mus_L3 = 0.1*8; % mm-1
mus_Bulk = 0.1*0.02; % mm-1

n=1.37;
g = 0.9;

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
if run_in_cluster == 1
    addpath('/pbs/home/c/ccaredda/private/mcx/utils');
    addpath('/pbs/home/c/ccaredda/private/mcxlab');
else
    addpath('/home/caredda/Soft/mcx/utils');
    addpath('/home/caredda/Soft/mcxlab');
    addpath('/home/caredda/Soft/iso2mesh');
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


% Set optical properties % [mua,mus,g,n]
% 0: Air
% 1: L2
% 2: L1
% 3: L3
% 4: Bulk
cfg.prop=[0 0 1 1; ...
mua_L1 mus_L1 g n ; ...
mua_L2 mus_L2 g n ; ...
mua_L3 mus_L3 g n ; ...
mua_Bulk mus_Bulk g n];


            
%Indexes of layers along z axis
N_1 = floor(L1_thickness_mm/model_resolution_in_mm);
N_2 = floor(L2_thickness_mm/model_resolution_in_mm);
N_3 = floor(L3_thickness_mm/model_resolution_in_mm);

% Create volume
cfg.vol = 4*ones(volume_square_size,volume_square_size,volume_square_size); %200 mm x 200 mm x 200 mm

%Configuration L2_L1_L3

% 0: Air
cfg.vol(:,:,1) = 0;
% L2
cfg.vol(:,:,2          :2+N_2-1) = 2; % L2
% 2: L1
cfg.vol(:,:,2+N_2      :2+N_1+N_2-1) = 1; % L1
% 3: L3
cfg.vol(:,:,2+N_1+N_2  :2+N_1+N_2+N_3-1) = 3; % L3
% 4: rest is L4

%Convert in uint8
cfg.vol = uint8(cfg.vol);

%Slice of the tissue
slice_tissue = squeeze(cfg.vol(floor(volume_square_size/2),floor(volume_square_size/2),:));

                        
%Init output
Diffuse_reflectance = zeros(8,Nb_measure);

%Process simulation with light source at the source
%position
cfg.srcpos=src_pos;


for n=1:Nb_measure
    % Random seed to obtain different results when running multiple simulations for the same input parameters
    cfg.seed = randi([0,99999],1);
    
    % calculate the fluence and partial path lengths
    [fluence,output_det]=mcxlab(cfg);
    
    % Detector pos
    for i=1:length(det_pos)
        %Diffuse reflectance
        Diffuse_reflectance(i,n) = fluence.dref(det_pos(i,1)+1,det_pos(i,2)+1,1);
    end
end


                        
%Save outputs
save('out_Phantom.mat','Diffuse_reflectance',...
                        'volume_square_size','src_pos',...
                        'det_pos','Lambdas','repetitions',...
                        'nphotons','model_resolution_in_mm',...
                        'L1_thickness_mm','L2_thickness_mm',...
                        'L3_thickness_mm','slice_tissue');

               
