
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
Lambda_array = [780, 810, 830, 840, 850, 890];

%Volume ize in mm
xdim_mm = 80;
ydim_mm = 80;
zdim_mm = 80;

%Source detector separation in mm
detectors_SD_mm = [30 35 40 45 50];

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


%Model fiber detector
cfg.radius_fiber_det_mm = 2.3;
cfg.reso_detector_mm = 0.1;


%Calculate optical properties for each layers
for p=1:length(SatO2_array)
    
    C_HbT = 35*1e-6;
    
    %Output
    Diffuse_reflectance = zeros(length(detectors_SD_mm), length(Lambda_array));
    Fluence_at_fiber_detector = zeros(length(detectors_SD_mm), length(Lambda_array)); 


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


        %Get fluence at detector pos
        %Init output
        [xi, yi, zi] = meshgrid(0.5:cfg.xdim_mm-0.5, 0.5:cfg.ydim_mm-0.5, 0.5:cfg.zdim_mm-0.5);
        Sensitivity_profile = zeros(length(xi), length(yi), length(zi), length(cfg.detectors_SD_mm));
    
    
        %Fluence at Fiber detector
        Phi_volume = griddata(cfg.node(:,1), cfg.node(:,2), cfg.node(:,3), phi(:,1), xi, yi, zi);
    
        %info interpolated image
        dim_interp_img = (2*ceil(cfg.radius_fiber_det_mm)+1)/cfg.reso_detector_mm;
    
            
    
        %Calculate sensitivity profile per detector
        for i=1:length(cfg.detectors_SD_mm)
    
            
            %Get fluence at fiber detector position
    
            %Detector is at position [ceil(radius_fiber_det_mm)+1 , ceil(radius_fiber_det_mm)+1]
            Phi_around_det = Phi_volume(cfg.detpos(i,2)+1 - ceil(cfg.radius_fiber_det_mm) : cfg.detpos(i,2)+1 + ceil(cfg.radius_fiber_det_mm) , ...
                                         cfg.detpos(i,1)+1 - ceil(cfg.radius_fiber_det_mm) : cfg.detpos(i,1)+1 + ceil(cfg.radius_fiber_det_mm));
            
            %Interpolate image around detector to match reso_detector_mm
            Rin  = imref2d(size(Phi_around_det), [0 2*ceil(cfg.radius_fiber_det_mm)], [0 2*ceil(cfg.radius_fiber_det_mm)]); 
            Rout = imref2d([dim_interp_img dim_interp_img],   [0 2*ceil(cfg.radius_fiber_det_mm)], [0 2*ceil(cfg.radius_fiber_det_mm)]);
            
            tform = affine2d(eye(3)); % identity (no geometric change, only resampling)
            J = imwarp(Phi_around_det, Rin, tform, ...
                       'OutputView', Rout, ...
                       'InterpolationMethod', 'cubic');
            
            % pixel center
            cx = (dim_interp_img+1)/2;  
            cy = (dim_interp_img+1)/2;
        
            %Mask on detector fiber position
            radius_fiber_det_px = cfg.radius_fiber_det_mm / Rout.PixelExtentInWorldX;
            [x, y] = meshgrid(1:dim_interp_img, 1:dim_interp_img);
            mask = (x - cx).^2 + (y - cy).^2 <= radius_fiber_det_px^2;
                
        
            Fluence_at_fiber_detector(i,i_Lambdas) = sum(J(mask));
    
        end

     
    end
    output_name = strcat(outdir,'/St_',num2str(SatO2_array(p)),'.mat');
    save(output_name,'Diffuse_reflectance','Fluence_at_fiber_detector');
end