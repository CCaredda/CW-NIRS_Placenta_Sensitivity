% TO DO calculate bananashape 
% get diffuse reflectance
clear
close all

run_in_cluster = 1;

% Repeat the simulation x times
repetitions = 1;
cfg.respin = repetitions;

% Number of photons
nphotons = 1e8;
cfg.nphoton = nphotons; 

% Maximum number of photons that can be detected 
cfg.maxdetphoton = 1e8;

%Size of 1 voxel in mm
model_resolution_in_mm = 1;

%Wavelength
Lambdas = 780;

%Arrays thickness
thickness_skin_in_mm = 3;
thickness_adipose_in_mm = 12;
thickness_muscle_in_mm = 20;



% Save diffuse reflectance
cfg.issaveref = 1;
                                           
%Acquisition time
cfg.tstart=0; % Starting time of the simulation (in seconds)
cfg.tend=1; % Ending time of the simulation (in seconds)
cfg.tstep=1; % Time-gate width of the simulation (in seconds)

% Calculate specular reflection if source is outside        
cfg.isspecular = 1; 
cfg.autopilot = 1;
   
% Voxel size in mm
cfg.unitinmm = model_resolution_in_mm;

% Add path
addpath('functions');
addpath('/home/caredda/Soft/mcx/utils');
addpath('/home/caredda/Soft/mcxlab');
addpath('/home/caredda/Soft/iso2mesh');



% Source from 0
cfg.issrcfrom0=1;

%Volume square size
volume_square_size = 200;

% Source type (values in accordance with VOlume size, se below)
cfg.srctype='pencil';
src_pos = [floor(volume_square_size/2 -1) floor((volume_square_size - 80/model_resolution_in_mm)/2 -1) volume_square_size-1];
cfg.srcdir=[0 0 -1];

cfg.srcpos=src_pos;

%Define multiple detectors  (values in accordance with VOlume size, se below)
det_pos = zeros(3,4);
offset_cm = 3;
for i=1:size(det_pos,1)
    det_pos(i,:) = [src_pos(1) src_pos(2)+(offset_cm+i)*floor(10/model_resolution_in_mm) volume_square_size 1]; % radius: 0.5 mm
end
cfg.detpos = det_pos;

cfg.savedetflag = 'dp'; %Save detector id and partial path length


% Set optical properties % [mua,mus,g,n]
% 0: Air
% 1: Skin
% 2: Adipose tissue
% 3: Muscle
% 4: Placenta
cfg.prop=[0 0 1 1; ...
0 0 1 1; ...
0 0 1 1 ; ...
0 0 1 1 ; ...
0 0 1 1];



%Indexes of layers along z axis
N_skin = floor(thickness_skin_in_mm/model_resolution_in_mm);
N_adipose = floor(thickness_adipose_in_mm/model_resolution_in_mm);
N_muscle = floor(thickness_muscle_in_mm/model_resolution_in_mm);

% Create volume
cfg.vol = 4*ones(volume_square_size,volume_square_size,volume_square_size); %200 mm x 200 mm x 200 mm

% 0: Air
cfg.vol(:,:,end) = 0;
% 1: skin
cfg.vol(:,:,end - N_skin : end-1) = 1; % Skin
% 2: Adipose tissue
cfg.vol(:,:,end - N_skin - N_adipose : end - N_skin -1) = 2; % Adipose
% 3: Muscle
cfg.vol(:,:,end - N_skin - N_adipose - N_muscle  :end - N_skin - N_adipose-1) = 3; % Muscle
% 4: rest is Placenta



mcxpreview(cfg)
axis off









