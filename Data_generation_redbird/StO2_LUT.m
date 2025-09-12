
clear
close all
clc

cluster = 1;

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
Lambda_array = [780, 810, 830, 840, 850, 890];

%Volume ize in mm
xdim_mm = 200;
ydim_mm = 200;
zdim_mm = 200;
% max_vol_mesh = [0.1; 0.1; 1; 1000];
max_vol_mesh = [0.5; 0.5; 1; 10000];

thickness_layers_mm_array = [1 2 7; ...
                             2 4 10; ...
                             3 5 12];

%Source detector separation in mm
detectors_SD_mm = [30 40 50];

%Saturation array
SatO2_muscle_array = [0.4, 0.4, 0.5, 0.7, 0.8, 0.8];
SatO2_placenta_array = [0.4, 0.8, 0.7, 0.5, 0.4, 0.8];


%Volume fraction of melanosome according to the color tones
% Modeling and Verification of Melanin Concentration on Human Skin Type
f_melanosome = [0.0255 0.155 0.305];


%Create out dir
outdir = strcat('StO2_LUT');
if ~exist(outdir, 'dir')
    mkdir(outdir)
end


subject=1;
    
clear cfg;

%Thickness layer
thickness_layers_mm = thickness_layers_mm_array(subject,:);

%Create 4 layers volume
cfg = create_meshed_volume_4layers(xdim_mm, ydim_mm, zdim_mm, thickness_layers_mm, max_vol_mesh, detectors_SD_mm, display);


%Calculate optical properties for each layers
for p=1:length(SatO2_muscle_array)
    for f=1:length(f_melanosome)
        C_HbT_muscle = 35*1e-6;
        C_HbT_placenta = 35*1e-6;
        SatO2_muscle = SatO2_muscle_array(p);
        SatO2_placenta = SatO2_placenta_array(p);
        f_mel = f_melanosome(f);

        %Output
        Diffuse_reflectance = zeros(length(detectors_SD_mm), length(Lambda_array));


        for i_Lambdas = 1:length(Lambda_array)
            % Calculate optical properties
            optical_prop = process_optical_properties_skin_Fat_muscle_placenta(Lambda_array(i_Lambdas),f_mel,SatO2_muscle, SatO2_placenta,C_HbT_muscle,C_HbT_placenta);


            %Get Sensitivity indexes for the 4 layers
            fprintf(1,strcat('Calculating sensitivity index\n'));
        
            %Calculate sensisitivity profile
            [sensitivity_profile, DR] = get_sensitivity_profiles(cfg, optical_prop);
            Diffuse_reflectance(:,i_Lambdas) = DR;
  
            output_name = strcat(outdir,'/out_subject_',num2str(subject),'St_muscle_',num2str(SatO2_muscle),'_St_placenta_',num2str(SatO2_placenta),'_fmel_',num2str(f_mel),'.mat');
            save(output_name,'Diffuse_reflectance');
        end
        
                            
    end
end