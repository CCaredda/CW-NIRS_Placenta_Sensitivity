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
% max_vol_mesh = [0.5; 0.5; 1; 10000];
max_vol_mesh = [1; 1; 1; 10000];

%Load subject info
thickness_skin_mm = ceil(readmatrix("../Subject_info/t_skin_cm.txt")*10);
thickness_adipose_mm = ceil(readmatrix("../Subject_info/t_adipose_cm.txt")*10);
thickness_muscle_mm = ceil(readmatrix("../Subject_info/t_adipose_cm.txt")*10);
Fitzpatrick_scale = ceil(readmatrix("../Subject_info/Fitzpatrick_scale.txt"));



%Source detector separation in mm
detectors_SD_mm = [30 40 50];

%Output
Sensitivity_indexes = zeros(length(detectors_SD_mm), 4, length(thickness_skin_mm));
Phi_detval = zeros(length(detectors_SD_mm), length(thickness_skin_mm));
DR_at_fiber_detector = zeros(length(detectors_SD_mm), length(thickness_skin_mm));
Fluence_at_fiber_detector = zeros(length(detectors_SD_mm), length(thickness_skin_mm));


for subject=1:length(thickness_skin_mm)
    fprintf(1,strcat('Subject',num2str(subject),'\n'));
    clear cfg;

    %Thickness layer
    thickness_layers_mm = [thickness_skin_mm(subject) thickness_adipose_mm(subject) thickness_muscle_mm(subject)];
    
    %Create 4 layers volume
    cfg = create_meshed_volume_4layers(xdim_mm, ydim_mm, zdim_mm, thickness_layers_mm, max_vol_mesh, detectors_SD_mm, display);

    %Model fiber detector
    cfg.radius_fiber_det_mm = 2.3;
    cfg.reso_detector_mm = 0.1;
   
    %Calculate optical properties for each layers
    C_HbT_muscle = 35*1e-6;
    C_HbT_placenta = 35*1e-6;
    SatO2_muscle = 0.6;
    SatO2_placenta = 0.8;
    f_mel = convert_Fitzpatrick_scale_f_mel(Fitzpatrick_scale(subject));

    fprintf(1,strcat('Tissue thickness ',num2str(thickness_layers_mm(1)),' ',num2str(thickness_layers_mm(2)),' ',num2str(thickness_layers_mm(3)),'\n'));
    fprintf(1,strcat('fmel ',num2str(f_mel),'\n'));


    % Calculate optical properties
    optical_prop = process_optical_properties_skin_Fat_muscle_placenta(Lambdas,f_mel,SatO2_muscle, SatO2_placenta,C_HbT_muscle,C_HbT_placenta);


    %Get Sensitivity indexes for the 4 layers
    fprintf(1,strcat('Calculating sensitivity index\n'));
    %Calculate sensisitivity profile
    [sensitivity_profile, Phi_det, Fluence_det, DR_det] = get_sensitivity_profiles(cfg, optical_prop);
    
    
    %Calculate sensiticity indexes
    S = get_sensitivity_index(cfg, sensitivity_profile, thickness_layers_mm);

    
    Sensitivity_indexes(:,:,subject) = S;
    Phi_detval(:,subject) = Phi_det;
    DR_at_fiber_detector(:,subject) = DR_det;
    Fluence_at_fiber_detector(:,subject) = Fluence_det;

    
    fprintf(1,strcat('End subject ',num2str(subject),'\n'));
end


save(strcat('Data_subjects_',num2str(Lambdas),'.mat'),'Phi_detval','Sensitivity_indexes', 'DR_at_fiber_detector', 'Fluence_at_fiber_detector');


