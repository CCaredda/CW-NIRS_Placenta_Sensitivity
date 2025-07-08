clear
close all
clc

tic

run_in_cluster = 0;

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
skin_thickness_array = [1,2,3,5];
adipose_thickness_array = [1,2,4,8,16];
muscle_thickness_array = [3,7,10,13,20,25];

%Saturation array
SatO2_muscle_array = 0.4:0.2:0.8;
SatO2_placenta_array = 0.4:0.2:0.8;


%HbT array
C_HbT_muscle_array = [15,25,35,50]*1e-6;
C_HbT_placenta_array = [15,25,35,50]*1e-6;


%Volume fraction of melanosome according to the color tones
% Modeling and Verification of Melanin Concentration on Human Skin Type
f_melanosome = [0.0255 0.155 0.305];





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
    % addpath('/home/caredda/Soft/mcx/utils');
    % addpath('/home/caredda/Soft/mcxlab');
    % addpath('/home/caredda/Soft/iso2mesh');
    addpath('/home/admin/Software/mcx-linux-x86_64-v2023/mcx/utils');
    addpath('/home/admin/Software/mcxlab-linux-x86_64-v2023/mcxlab');  
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
if ~exist('output_lookup_table', 'dir')
    mkdir('output_lookup_table')
end
save('output_lookup_table/cst.mat','volume_square_size','src_pos','det_pos','Lambdas','repetitions','nphotons','model_resolution_in_mm','SatO2_placenta_array','SatO2_muscle_array','skin_thickness_array','adipose_thickness_array','muscle_thickness_array','f_melanosome');


%Loop on HbT muscle
for C_HbT_muscle = C_HbT_muscle_array
    for C_HbT_placenta = C_HbT_placenta_array
        % Loop on SatO2 muscle
        for SatO2_muscle = SatO2_muscle_array
            % Loop on SatO2 placenta
            for SatO2_placenta = SatO2_placenta_array
                %Loop on color tones
                for f_mel = f_melanosome
                    
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
        
        
                    %Loop over skin thickness
                    for thickness_skin_in_mm = skin_thickness_array
                       
                        %Loop over Fat thickness
                        for thickness_adipose_in_mm = adipose_thickness_array
        
                            %Loop over muscle thickness
                            for thickness_muscle_in_mm = muscle_thickness_array
            
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
                                Diffuse_reflectance = zeros(8,1);
                                Sensitivity_indexes = zeros(8,4);
                                Nb_detected_photons = zeros(8,1);
                                
        
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
                                    %Nb of detected photons
                                    Nb_detected_photons(i) = length(find(output_det.detid==i));
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
        
        
                                    %Get sensitivity indexes
                                    for m=1:4
                                        % Find indices of tissue m
                                        indices = find(slice_tissue == m);
                                        % Compute sensitivity
                                        Sensitivity_indexes(d,m) = sum(Sensitivity_profile(:, :, indices),"all");
                                    end
                                end
                            
                                %Save outputs
                                output_name = strcat('output_lookup_table/out_St_muscle_',num2str(SatO2_muscle),'_St_placenta_',num2str(SatO2_placenta),'_Thick_skin_',num2str(thickness_skin_in_mm),'_Thick_adipose_',num2str(thickness_adipose_in_mm),'_Thick_muscle_',num2str(thickness_muscle_in_mm),'f_mel',num2str(f_mel),'_HbT_muscle_umol_',num2str(C_HbT_muscle*1e6),'_HbT_placenta_umol_',num2str(C_HbT_placenta*1e6),'.mat');
                                save(output_name,'Diffuse_reflectance','Sensitivity_indexes',"Nb_detected_photons");
        
                         
                            end
                        end
                    end
                end
            end
        end
    end 
end


toc








