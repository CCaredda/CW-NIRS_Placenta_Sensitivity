
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

%Source detector separation in mm
detectors_SD_mm = [30 40 50];

%Saturation array
SatO2_array = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];


%Create out dir
outdir = strcat('StO2_semi_ininite');
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

clear cfg;

%Create homogeneous volume
[cfg.node, cfg.face, cfg.elem] = meshabox([0 0 0], [xdim_mm ydim_mm zdim_mm], 1);
cfg.seg = ones(size(cfg.elem, 1), 1);

%Light source
src_pos = [(xdim_mm/2-1) (ydim_mm - detectors_SD_mm(end))/2-1 0];
cfg.srcdir=[0 0 1];
cfg.srcpos = src_pos;
cfg.srctype = 'pencil';

%Detectors
%Do not model the fiber optics (it is scaled with the calibration procedure)
cfg.detdir = [0 0 1];
det_pos = zeros(length(detectors_SD_mm),3);
for i=1:length(detectors_SD_mm)
    det_pos(i,:) = [src_pos(1) src_pos(2)+detectors_SD_mm(i) 0]; % No radius!!
end
cfg.detpos = det_pos;


%Calculate optical properties for each layers
for p=1:length(SatO2_array)
    
    C_HbT = 35*1e-6;
    
    %Output
    Diffuse_reflectance = zeros(length(detectors_SD_mm), length(Lambda_array));


    for i_Lambdas = 1:length(Lambda_array)
        % Calculate optical properties
        optical_prop = process_optical_properties_skin_Fat_muscle_placenta(Lambda_array(i_Lambdas),0,0, SatO2_array(p),0,C_HbT);

        cfg.prop=[0 0 1 1; ...
                  optical_prop.mua_placenta optical_prop.mu_s_placenta optical_prop.g_placenta optical_prop.n_placenta];


        %Prepare all missing values in cfg struct
        % disp('Forward solver...')
        cfg = rbmeshprep(cfg);
    
    
        %Forward solver
        [DR, phi] = rbrun(cfg);
        Diffuse_reflectance(:,i_Lambdas) = DR;


        %Process srs
        h = 6.3e-4;
srs_mua_30_50_c1(lambda) = 1 / (3 * (1 - (h * lambda_all(lambda)))) *  (log(10) * (srs_ss_30_50_c1(lambda))-(2/mean([30 50]))).^2;
srs_mua_30_50_c2(lambda) = 1 / (3 * (1 - (h * lambda_all(lambda)))) *  (log(10) * (srs_ss_30_50_c2(lambda))-(2/mean([30 50]))).^2;


        output_name = strcat(outdir,'/St_',num2str(SatO2_array(p)),'.mat');
        save(output_name,'Diffuse_reflectance');
    end
end