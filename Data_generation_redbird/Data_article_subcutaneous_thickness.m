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
% Lambdas_array = [780 840 890];
Lambdas_array = [780];


%Volume ize in mm
xdim_mm = 200;
ydim_mm = 200;
zdim_mm = 200;
% max_vol_mesh = [0.1; 0.1; 1; 1000];
% max_vol_mesh = [0.5; 0.5; 1; 10000];
max_vol_mesh = [1; 3; 3; 10000];

% thickness_layers_mm_array = [2 2 7; ...
%                              2 4 10; ...
%                              2 5 12; ...
%                              2 7 17];

thickness_layers_mm_array = [1 5 12; ...
                             2 5 12; ...
                             3 5 12];





%Source detector separation in mm
detectors_SD_mm = [30 40 50];



%Create out dir
outdir = 'data_article_subcutaneous';
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
    
    
for subject=1:size(thickness_layers_mm_array,1)
    
    clear cfg;

    %Thickness layer
    thickness_layers_mm = thickness_layers_mm_array(subject,:);

    %Create 4 layers volume
    cfg = create_meshed_volume_4layers(xdim_mm, ydim_mm, zdim_mm, thickness_layers_mm, max_vol_mesh, detectors_SD_mm, display);

    %Model fiber detector
    cfg.radius_fiber_det_mm = 2.3;
    cfg.reso_detector_mm = 0.1;


    for Lambdas = Lambdas_array
        disp(Lambdas)

        %Calculate optical properties for each layers

        C_HbT_muscle = 35*1e-6;
        C_HbT_placenta = 35*1e-6;
        SatO2_muscle = 0.6;
        SatO2_placenta = 0.8;
        f_mel = 0.0255;
        % f_mel = 0.305;

        % Calculate optical properties
        optical_prop = process_optical_properties_skin_Fat_muscle_placenta(Lambdas,f_mel,SatO2_muscle, SatO2_placenta,C_HbT_muscle,C_HbT_placenta);


        %Get Sensitivity indexes for the 4 layers
        fprintf(1,strcat('Calculating sensitivity index\n'));
        
        %Calculate sensisitivity profile
        [sensitivity_profile, Diffuse_reflectance, Fluence_at_fiber_detector] = get_sensitivity_profiles(cfg, optical_prop);

        
        %Calculate sensiticity indexes
        Sensitivity_indexes = get_sensitivity_index(cfg, sensitivity_profile, thickness_layers_mm);

        output_name = strcat(outdir,'/out_',num2str(Lambdas),'_St_muscle_',num2str(SatO2_muscle),'_St_placenta_',num2str(SatO2_placenta),'_Thick_skin_',num2str(thickness_layers_mm(1)),'_Thick_adipose_',num2str(thickness_layers_mm(2)),'_Thick_muscle_',num2str(thickness_layers_mm(3)),'f_mel',num2str(f_mel),'_HbT_muscle_umol_',num2str(C_HbT_muscle*1e6),'_HbT_placenta_umol_',num2str(C_HbT_placenta*1e6),'.mat');
        save(output_name,'Diffuse_reflectance','Sensitivity_indexes', 'Fluence_at_fiber_detector');
                            
    end
end

