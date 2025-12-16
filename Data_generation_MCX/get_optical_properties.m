%% MCX


clear
close all
clc

addpath("../functions/");

% Fitzpatrick_scale = ceil(readmatrix("../Subject_info/Fitzpatrick_scale.txt"));



%Saturation array
C_HbT_placenta_array = [15,25,35,50]*1e-6;

%Volume fraction of melanosome according to the color tones
% Modeling and Verification of Melanin Concentration on Human Skin Type
f_melanosome = [0.0255 0.155 0.305];



%Create out dir
outdir = 'optical_properties';
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

for Lambdas = [780,840,890]

    %Init output
    mua_skin = zeros(length(C_HbT_placenta_array)*length(f_melanosome),1);
    mu_s_skin = zeros(size(mua_skin));
    mua_adipose = zeros(size(mua_skin)); 
    mu_s_adipose = zeros(size(mua_skin));
    mua_muscle = zeros(size(mua_skin));
    mu_s_muscle = zeros(size(mua_skin));
    mua_placenta = zeros(size(mua_skin));
    mu_s_placenta = zeros(size(mua_skin));


    %Calculate optical properties for each layers
    id = 1;
    for p=1:length(C_HbT_placenta_array)
        for f=1:length(f_melanosome)
            C_HbT_muscle = 35*1e-6;
            C_HbT_placenta = C_HbT_placenta_array(p);
            SatO2_muscle = 0.6;
            SatO2_placenta = 0.8;
            f_mel = f_melanosome(f);
            
    
            % Define optical properties
            optical_prop = process_optical_properties_skin_Fat_muscle_placenta(Lambdas,f_mel,SatO2_muscle, SatO2_placenta,C_HbT_muscle,C_HbT_placenta);
    
            mua_skin(id) = optical_prop.mua_skin;
            mu_s_skin(id) = optical_prop.mu_s_skin;
            mua_adipose(id) = optical_prop.mua_adipose;
            mu_s_adipose(id) = optical_prop.mu_s_adipose;
            mua_muscle(id) = optical_prop.mua_muscle;
            mu_s_muscle(id) = optical_prop.mu_s_muscle;
            mua_placenta(id) = optical_prop.mua_placenta;
            mu_s_placenta(id) = optical_prop.mu_s_placenta;
    
            id = id+1;
    
        end
    end
    
    %Save outputs
    output_name = strcat(outdir,'/optical_prop_data_',num2str(Lambdas),'.mat');
    save(output_name,'C_HbT_placenta_array','f_melanosome','mua_skin','mu_s_skin','mua_adipose','mu_s_adipose','mua_muscle','mu_s_muscle','mua_placenta','mu_s_placenta');
         


    % %Data subjects
    % 
    % %Init output
    % mua_skin = zeros(length(Fitzpatrick_scale),1);
    % mu_s_skin = zeros(size(mua_skin));
    % mua_adipose = zeros(size(mua_skin)); 
    % mu_s_adipose = zeros(size(mua_skin));
    % mua_muscle = zeros(size(mua_skin));
    % mu_s_muscle = zeros(size(mua_skin));
    % mua_placenta = zeros(size(mua_skin));
    % mu_s_placenta = zeros(size(mua_skin));
    % 
    % for subject=1:length(Fitzpatrick_scale)
    %     %Calculate optical properties for each layers
    %     C_HbT_muscle = 35*1e-6;
    %     C_HbT_placenta = 35*1e-6;
    %     SatO2_muscle = 0.6;
    %     SatO2_placenta = 0.8;
    %     f_mel = convert_Fitzpatrick_scale_f_mel(Fitzpatrick_scale(subject));
    % 
    %     % fprintf(1,strcat('Tissue thickness ',num2str(thickness_layers_mm(1)),' ',num2str(thickness_layers_mm(2)),' ',num2str(thickness_layers_mm(3)),'\n'));
    %     % fprintf(1,strcat('fmel ',num2str(f_mel),'\n'));
    % 
    % 
    %     % Calculate optical properties
    %     optical_prop = process_optical_properties_skin_Fat_muscle_placenta(Lambdas,f_mel,SatO2_muscle, SatO2_placenta,C_HbT_muscle,C_HbT_placenta);
    % 
    %     mua_skin(subject) = optical_prop.mua_skin;
    %     mu_s_skin(subject) = optical_prop.mu_s_skin;
    %     mua_adipose(subject) = optical_prop.mua_adipose;
    %     mu_s_adipose(subject) = optical_prop.mu_s_adipose;
    %     mua_muscle(subject) = optical_prop.mua_muscle;
    %     mu_s_muscle(subject) = optical_prop.mu_s_muscle;
    %     mua_placenta(subject) = optical_prop.mua_placenta;
    %     mu_s_placenta(subject) = optical_prop.mu_s_placenta;
    % 
    % end
    % 
    % %Save outputs
    % output_name = strcat(outdir,'/optical_prop_subjects_',num2str(Lambdas),'.mat');
    % save(output_name,'C_HbT_placenta_array','f_melanosome','mua_skin','mu_s_skin','mua_adipose','mu_s_adipose','mua_muscle','mu_s_muscle','mua_placenta','mu_s_placenta');
    % 


end




%%


