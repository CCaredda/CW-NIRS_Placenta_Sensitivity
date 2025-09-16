clear
close all
clc

addpath('~/Soft/redbird/matlab');
addpath('~/Soft/iso2mesh');
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


%Fiber detector size
radius_fiber_det_mm = 2.3;
reso_detector_mm = 0.1;


%Thickness layer
thickness_layers_mm = [2 4 10];

%Source detector separation in mm
detectors_SD_mm = [30 40 50];

clear cfg;

%Create 4 layers volume
cfg = create_meshed_volume_4layers(xdim_mm, ydim_mm, zdim_mm, thickness_layers_mm, max_vol_mesh, detectors_SD_mm, display);

C_HbT_muscle = 25*1e-6;
C_HbT_placenta = 25*1e-6;
SatO2_muscle = 0.6;
SatO2_placenta = 0.8;
% f_mel = 0.0255;
f_mel = 0.305;

% Calculate optical properties
optical_prop = process_optical_properties_skin_Fat_muscle_placenta(Lambdas,f_mel,SatO2_muscle, SatO2_placenta,C_HbT_muscle,C_HbT_placenta);


cfg.prop=[0 0 1 1; ...
optical_prop.mua_skin optical_prop.mu_s_skin optical_prop.g_skin optical_prop.n_skin ; ...
optical_prop.mua_adipose optical_prop.mu_s_adipose optical_prop.g_adipose optical_prop.n_adipose ; ...
optical_prop.mua_muscle optical_prop.mu_s_muscle optical_prop.g_muscle optical_prop.n_muscle ; ...
optical_prop.mua_placenta optical_prop.mu_s_placenta optical_prop.g_placenta optical_prop.n_placenta];


%Prepare all missing values in cfg struct
% disp('Forward solver...')
cfg = rbmeshprep(cfg);


%Forward solver
[Diffuse_reflectance, phi] = rbrun(cfg);
%Set negative values to 0
phi = abs(phi);

%interpolate volume
[xi, yi, zi] = meshgrid(0.5:xdim_mm-0.5, 0.5:ydim_mm-0.5, 0.5:zdim_mm-0.5);
Phi_volume = griddata(cfg.node(:,1), cfg.node(:,2), cfg.node(:,3), phi(:,1), xi, yi, zi);



%info interpolated image
dim_interp_img = (2*ceil(radius_fiber_det_mm)+1)/reso_detector_mm;




DR_phi = [];
DR_phi_test = [];



%Print detvalue at detpos
for i=1:length(detectors_SD_mm)
    %Store fluence values at detector
    DR_phi(end+1) = Phi_volume(cfg.detpos(i,2)+1,cfg.detpos(i,1)+1);

    %Select image around detector
    %Detector is at position [ceil(radius_fiber_det_mm)+1 , ceil(radius_fiber_det_mm)+1]
    Phi_around_det = Phi_volume(cfg.detpos(i,2)+1 - ceil(radius_fiber_det_mm) : cfg.detpos(i,2)+1 + ceil(radius_fiber_det_mm) , ...
                                 cfg.detpos(i,1)+1 - ceil(radius_fiber_det_mm) : cfg.detpos(i,1)+1 + ceil(radius_fiber_det_mm));
    
    %Interpolate image around detector to match reso_detector_mm
    Rin  = imref2d(size(Phi_around_det), [0 2*ceil(radius_fiber_det_mm)], [0 2*ceil(radius_fiber_det_mm)]); 
    Rout = imref2d([dim_interp_img dim_interp_img],   [0 2*ceil(radius_fiber_det_mm)], [0 2*ceil(radius_fiber_det_mm)]);
    
    tform = affine2d(eye(3)); % identity (no geometric change, only resampling)
    J = imwarp(Phi_around_det, Rin, tform, ...
               'OutputView', Rout, ...
               'InterpolationMethod', 'cubic');
    
    % pixel center
    cx = (dim_interp_img+1)/2;  
    cy = (dim_interp_img+1)/2;

    %Mask on detector fiber position
    radius_fiber_det_px = radius_fiber_det_mm / Rout.PixelExtentInWorldX;
    [x, y] = meshgrid(1:dim_interp_img, 1:dim_interp_img);
    mask = (x - cx).^2 + (y - cy).^2 <= radius_fiber_det_px^2;
        
    % DR_phi_test(end+1) = J(floor(cx),floor(cy));

    DR_phi_test(end+1) = sum(J(mask));
    % DR_phi_test(end+1) = mean(J(mask));
end


close('all')
figure(1)
hold on
plot(Diffuse_reflectance,'r')
plot(DR_phi,'b')
plot(DR_phi_test,'g')
