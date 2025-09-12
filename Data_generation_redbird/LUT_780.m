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



%Wavelength
Lambdas = 780;
% Lambdas = 840;
% Lambdas = 890;


%Arrays thickness
skin_thickness_array_mm = [1,2,3,4,5];
adipose_thickness_array_mm = [1,2,4,8,16];
muscle_thickness_array_mm = [3,7,10,13,20,25];

%Saturation array
SatO2_muscle_array = 0.4:0.2:0.8;
SatO2_placenta_array = 0.4:0.2:0.8;


%HbT array
C_HbT_muscle_array = [15,25,35,50]*1e-6;
C_HbT_placenta_array = [15,25,35,50]*1e-6;


%Volume fraction of melanosome according to the color tones
% Modeling and Verification of Melanin Concentration on Human Skin Type
f_melanosome = [0.0255 0.155 0.305];


%Volume ize in mm
xdim_mm = 200;
ydim_mm = 200;
zdim_mm = 200;
% max_vol_mesh = [0.1; 0.1; 1; 1000];
max_vol_mesh = [0.5; 0.5; 1; 10000];


%Info detectors
detectors_SD_mm = [30, 40, 50];

%Save constant
outdir = strcat('output_lookup_table_',num2str(Lambdas));
if ~exist(outdir, 'dir')
    mkdir(outdir)
end


%Define Volume
for thickness_skin_in_mm = skin_thickness_array_mm
    fprintf(1,strcat('Skin thickness: ',num2str(thickness_skin_in_mm),'\n'));
    
    for thickness_adipose_in_mm = adipose_thickness_array_mm
        fprintf(1,strcat('Adipose thickness: ',num2str(thickness_adipose_in_mm),'\n'));

        for thickness_muscle_in_mm = muscle_thickness_array_mm
            fprintf(1,strcat('Muscle thickness: ',num2str(thickness_muscle_in_mm),'\n'));

            %Thickness layer
            thickness_layers_mm = [thickness_skin_in_mm thickness_adipose_in_mm thickness_muscle_in_mm];
    
            %Create 4 layers volume
            cfg = create_meshed_volume_4layers(xdim_mm, ydim_mm, zdim_mm, thickness_layers_mm, max_vol_mesh, detectors_SD_mm, display);
        
                    
            %Define optical properties for all wavelengths
            for C_HbT_muscle = C_HbT_muscle_array
                for C_HbT_placenta = C_HbT_placenta_array
                    for SatO2_muscle = SatO2_muscle_array
                        for SatO2_placenta = SatO2_placenta_array
                            for f_mel = f_melanosome
            

                                % Calculate optical properties
                                optical_prop = process_optical_properties_skin_Fat_muscle_placenta(Lambdas,f_mel,SatO2_muscle, SatO2_placenta,C_HbT_muscle,C_HbT_placenta);
                                
                                %Calculate sensisitivity profile
                                [sensitivity_profile, Diffuse_reflectance] = get_sensitivity_profiles(cfg, optical_prop);
                                
                                %Calculate sensiticity indexes
                                Sensitivity_indexes = get_sensitivity_index(cfg, sensitivity_profile, thickness_layers_mm);

                                %Save outputs
                                output_name = strcat(outdir,'/out_St_muscle_',num2str(SatO2_muscle),'_St_placenta_',num2str(SatO2_placenta),'_Thick_skin_',num2str(thickness_skin_in_mm),'_Thick_adipose_',num2str(thickness_adipose_in_mm),'_Thick_muscle_',num2str(thickness_muscle_in_mm),'f_mel',num2str(f_mel),'_HbT_muscle_umol_',num2str(C_HbT_muscle*1e6),'_HbT_placenta_umol_',num2str(C_HbT_placenta*1e6),'.mat');
                                save(output_name,'Diffuse_reflectance','Sensitivity_indexes');
                                
                            end
                        end
                    end
                end
            end
        end
    end
end

            
            
            






















