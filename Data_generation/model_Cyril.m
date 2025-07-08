% To do: plot bananashap maps
% Check scattering and absorption coeff with Fred code


clear
close all

SatO2_muscle = 0.6;
SatO2_placenta = 0.4;

model_resolution_in_mm = 2;

Lambdas = 710:10:900;
% Lambdas = 780;

run_in_cluster = 0;
nb_repeat = 1; %Nb of repetitions used in MCX
nb_photons = 1e9;%1e6;

out_path = 'output/';

% Add path
if run_in_cluster == 1
addpath('functions');
    addpath('/pbs/home/c/ccaredda/private/mcx/utils');
    addpath('/pbs/home/c/ccaredda/private/mcxlab');
else
    addpath('/home/caredda/Soft/mcx/utils');
    addpath('/home/caredda/Soft/mcxlab');
    % addpath('/home/admin/Software/mcx-linux-x86_64-v2023/mcx/utils');
    % addpath('/home/admin/Software/mcxlab-linux-x86_64-v2023/mcxlab');  
end

% Number of photons
cfg.nphoton=nb_photons; 

% Repeat the simulation x times
cfg.respin=nb_repeat;

% Maximum number of photons that can be detected 
cfg.maxdetphoton = 1e8;

% GPU processing
cfg.gpuid=1; 

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

% Create volume
% cfg.vol = 3*ones(60,60,60); %120 mm x 120 mm x 120 mm
cfg.vol = 3*ones(60,120,60); %120 mm x 120 mm x 120 mm

% 0: Air
cfg.vol(:,:,1) = 0;
% 1: Adipose tissue
cfg.vol(:,:,2:3) = 1; % 4mm adipose
% 2: Muscle
cfg.vol(:,:,4:9) = 2; % 10mm muscle
% 3: rest is Placenta

%Convert in uint8
cfg.vol = uint8(cfg.vol);

slice_tissue = cfg.vol(39,60,:);
save('slice_tissue.mat','slice_tissue');

% Source from 0
cfg.issrcfrom0=1;

% Source type 
cfg.srctype='pencil';
% cfg.srcpos=[29 9 0];
cfg.srcpos=[29 39 0];
cfg.srcdir=[0 0 1];

% Boundary conditions
% cfg.bc='ccrccr';
% ccrccc: Cyclic BC except for the top and bottom face (Fresnel reflection)
    
%Define multiple detectors
det_pos = zeros(8,4);
for i=1:8
    % det_pos(i,:) = [29 14+(i-1)*5 1 0.5]; % radius: 0.5 mm
    det_pos(i,:) = [29 44+(i-1)*5 1 0.5]; % radius: 0.5 mm
end
cfg.detpos = det_pos;

cfg.savedetflag = 'dp'; %Save detector id and partial path length

% Define optical properties
[mua_adipose,mua_muscle,mua_placenta,mu_s_adipose,mu_s_muscle,mu_s_placenta,g,n] = process_optical_properties(Lambdas, SatO2_muscle, SatO2_placenta);

%Init output
Diffuse_reflectance = zeros(8,length(Lambdas));
Mean_path = zeros(8,length(Lambdas));

% process simulations
for l=1:length(Lambdas)
    disp(strcat("Simulation lambda ",num2str(Lambdas(l))))
        
    % Set optical properties % [mua,mus,g,n]
    % 0: Air
    % 1: Adipose tissue
    % 2: Muscle
    % 3: Placenta

    cfg.prop=[0 0 1 1; ...
    mua_adipose(l) mu_s_adipose(l) g n ; ...
    mua_muscle(l) mu_s_muscle(l) g n ; ...
    mua_placenta(l) mu_s_placenta(l) g n];

    % Random seed to obtain different results when running multiple simulations for the same input parameters
    cfg.seed = randi([0,99999],1);
    
    % calculate the fluence and partial path lengths
    % [flux,output_det]=mcxlab(cfg);
    [fluence,output_det,vol,seed,trajectory]=mcxlab(cfg);
      
    % Detector pos
    for i=1:length(det_pos)
        %Diffuse reflectance
        Diffuse_reflectance(i,l) = fluence.dref(det_pos(i,1)+1,det_pos(i,2)+1,1);
        
        %Mean path
        Mean_path(i,l) = sum(mcxmeanpath_detID(i,output_det));
    
    end

     
end

% Save results

save(strcat('out_muscle_',num2str(SatO2_muscle),'_placenta_',num2str(SatO2_placenta),'.mat'),'Diffuse_reflectance','Mean_path','Lambdas','slice_tissue');

