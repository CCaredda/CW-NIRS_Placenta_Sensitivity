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

thickness_layers_mm_array = [1 2 7; ...
                             2 4 10; ...
                             3 5 12];

%Source detector separation in mm
detectors_SD_mm = [30 40 50];

%Saturation array
C_HbT_placenta_array = [15,25,35,50]*1e-6;

%Volume fraction of melanosome according to the color tones
% Modeling and Verification of Melanin Concentration on Human Skin Type
f_melanosome = [0.0255 0.155 0.305];


%Create out dir
outdir = strcat('data_article_',num2str(Lambdas));
if ~exist(outdir, 'dir')
    mkdir(outdir)
end


for subject=1:size(thickness_layers_mm_array,1)
    
    clear cfg;

    %Thickness layer
    thickness_layers_mm = thickness_layers_mm_array(subject,:);

    %Create 4 layers volume
    cfg = create_meshed_volume_4layers(xdim_mm, ydim_mm, zdim_mm, thickness_layers_mm, max_vol_mesh, detectors_SD_mm, display);


    %Calculate optical properties for each layers
    for p=1:length(C_HbT_placenta_array)
        for f=1:length(f_melanosome)
            C_HbT_muscle = 25*1e-6;
            C_HbT_placenta = C_HbT_placenta_array(p);
            SatO2_muscle = 0.6;
            SatO2_placenta = 0.8;
            f_mel = f_melanosome(f);

            % Calculate optical properties
            optical_prop = process_optical_properties_skin_Fat_muscle_placenta(Lambdas,f_mel,SatO2_muscle, SatO2_placenta,C_HbT_muscle,C_HbT_placenta);


            %Get Sensitivity indexes for the 4 layers
            fprintf(1,strcat('Calculating sensitivity index\n'));
            
            %Calculate sensisitivity profile
            [sensitivity_profile, Diffuse_reflectance] = get_sensitivity_profiles(cfg, optical_prop);
            
            %Calculate sensiticity indexes
            Sensitivity_indexes = get_sensitivity_index(cfg, sensitivity_profile, thickness_layers_mm);

            output_name = strcat(outdir,'/out_St_muscle_',num2str(SatO2_muscle),'_St_placenta_',num2str(SatO2_placenta),'_Thick_skin_',num2str(thickness_layers_mm(1)),'_Thick_adipose_',num2str(thickness_layers_mm(2)),'_Thick_muscle_',num2str(thickness_layers_mm(3)),'f_mel',num2str(f_mel),'_HbT_muscle_umol_',num2str(C_HbT_muscle*1e6),'_HbT_placenta_umol_',num2str(C_HbT_placenta*1e6),'.mat');
            save(output_name,'Diffuse_reflectance','Sensitivity_indexes');
                                
        end
    end
end
