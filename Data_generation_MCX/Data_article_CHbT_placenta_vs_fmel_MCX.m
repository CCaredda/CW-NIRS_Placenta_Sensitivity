%% MCX


clear
close all
clc

subject = 1;

% Number of photons
nphotons = 1e9;
% nphotons = 1e2;


%Wavelength
Lambdas = 780;

%Create out dir
outdir = strcat('data_article_MCX_',num2str(Lambdas));
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

%Thickness layers (in mm)
thickness_layers_mm_array = [1 2 7; ...
                             2 4 10; ...
                             3 5 12];

%Model fiber detector
radius_fiber_det_mm = 2.3;

%Source detector separation in mm
detectors_SD_mm = [30 40 50];

%Saturation array
C_HbT_placenta_array = [15,25,35,50]*1e-6;

%Volume fraction of melanosome according to the color tones
% Modeling and Verification of Melanin Concentration on Human Skin Type
f_melanosome = [0.0255 0.155 0.305];

run_in_cluster = 0;

% Repeat the simulation x times
repetitions = 1;
cfg.respin = repetitions;

% Number of photons
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

% Add path
addpath('../functions');
if run_in_cluster == 1
    addpath('~/private/mcx/utils');
    addpath('~/private/mcxlab');
else
    addpath('~/Soft/mcx/utils');
    addpath('~/Soft/mcxlab');
    addpath('~/Soft/iso2mesh'); 
end



% Source from 0
cfg.issrcfrom0=1;

%Volume square size
volume_square_size = 200;

%Volume ize in mm
xdim_mm = volume_square_size;
ydim_mm = volume_square_size;
zdim_mm = volume_square_size;

%Light source
src_pos = [(xdim_mm/2-1) (ydim_mm - detectors_SD_mm(end))/2-1 0];
cfg.srcdir=[0 0 1];
cfg.srctype = 'pencil';

%Detectors
det_pos = zeros(length(detectors_SD_mm),4);
for i=1:length(detectors_SD_mm)
    det_pos(i,:) = [src_pos(1) src_pos(2)+detectors_SD_mm(i) 1 radius_fiber_det_mm];
end
cfg.detpos = det_pos;


cfg.savedetflag = 'dp'; %Save detector id and partial path length

% GPU processing
cfg.gpuid=1;


%Thickness layer
thickness_layers_mm = thickness_layers_mm_array(subject,:);

%Calculate optical properties for each layers
for p=1:length(C_HbT_placenta_array)
    for f=1:length(f_melanosome)
        C_HbT_muscle = 25*1e-6;
        C_HbT_placenta = C_HbT_placenta_array(p);
        SatO2_muscle = 0.6;
        SatO2_placenta = 0.8;
        f_mel = f_melanosome(f);
        
        %Tissue thickness
        thickness_skin_in_mm = thickness_layers_mm(1);
        thickness_adipose_in_mm = thickness_layers_mm(2);
        thickness_muscle_in_mm = thickness_layers_mm(3);

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
        N_skin = floor(thickness_layers_mm(1)/model_resolution_in_mm);
        N_adipose = floor(thickness_layers_mm(2)/model_resolution_in_mm);
        N_muscle = floor(thickness_layers_mm(3)/model_resolution_in_mm);
        
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
        Diffuse_reflectance = zeros(length(detectors_SD_mm),1);
        Sensitivity_indexes = zeros(length(detectors_SD_mm),4);
        

        %Process simulation with light source at the source
        %position
        cfg.srcpos=src_pos;

        % Random seed to obtain different results when running multiple simulations for the same input parameters
        cfg.seed = randi([0,99999],1);
        
        % calculate the fluence and partial path lengths
        [fluence,output_det]=mcxlab(cfg);
          
        % Detector pos
        for i=1:length(detectors_SD_mm)
            %Diffuse reflectance
            Diffuse_reflectance(i) = fluence.dref(det_pos(i,1)+1,det_pos(i,2)+1,1);
        end
        
        
        
        % Calculate diffuse reflectance at fiber detector location using fluence
        [DR_at_fiber_detector,Fluence_at_fiber_detector] = get_Diffuse_reflectance_from_Fluence_MCX(fluence.data, det_pos, radius_fiber_det_mm, 0.1, optical_prop);

        %Add numbers of detected photons

        
        % Compute simulation at detector location
        for d = 1:length(detectors_SD_mm)
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
        output_name = strcat(outdir,'/out_St_muscle_',num2str(SatO2_muscle),'_St_placenta_',num2str(SatO2_placenta),'_Thick_skin_',num2str(thickness_skin_in_mm),'_Thick_adipose_',num2str(thickness_adipose_in_mm),'_Thick_muscle_',num2str(thickness_muscle_in_mm),'f_mel',num2str(f_mel),'_HbT_muscle_umol_',num2str(C_HbT_muscle*1e6),'_HbT_placenta_umol_',num2str(C_HbT_placenta*1e6),'.mat');
        save(output_name,'Diffuse_reflectance','Sensitivity_indexes', 'DR_at_fiber_detector', 'Fluence_at_fiber_detector');

                 
    end
end
