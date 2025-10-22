function [Sensitivity_profile, phi_detval, DR_at_fiber_detector] = get_sensitivity_profiles(cfg, optical_prop)
    cfg.prop=[0 0 1 1; ...
    optical_prop.mua_skin optical_prop.mu_s_skin optical_prop.g_skin optical_prop.n_skin ; ...
    optical_prop.mua_adipose optical_prop.mu_s_adipose optical_prop.g_adipose optical_prop.n_adipose ; ...
    optical_prop.mua_muscle optical_prop.mu_s_muscle optical_prop.g_muscle optical_prop.n_muscle ; ...
    optical_prop.mua_placenta optical_prop.mu_s_placenta optical_prop.g_placenta optical_prop.n_placenta];
    

    %Prepare all missing values in cfg struct
    % disp('Forward solver...')
    cfg = rbmeshprep(cfg);
    
    
    %Forward solver
    [phi_detval, phi] = rbrun(cfg);
    
    
    
    %Set negative values to 0
    phi = abs(phi);


    %Init output
    [xi, yi, zi] = meshgrid(0.5:cfg.xdim_mm-0.5, 0.5:cfg.ydim_mm-0.5, 0.5:cfg.zdim_mm-0.5);
    Sensitivity_profile = zeros(length(xi), length(yi), length(zi), length(cfg.detectors_SD_mm));


    %Fluence at Fiber detector
    DR_at_fiber_detector = [];
    Phi_volume = griddata(cfg.node(:,1), cfg.node(:,2), cfg.node(:,3), phi(:,1), xi, yi, zi);
    
    %Calculate derivative of fluence along z axis (resolution is 1mm3)
    [Fx,Fy,Fz] = gradient(Phi_volume,1,1,1);
    
    %info interpolated image
    dim_interp_img = (2*ceil(cfg.radius_fiber_det_mm)+1)/cfg.reso_detector_mm;


    %Calculate sensitivity profile per detector
    for i=1:length(cfg.detectors_SD_mm)

        %Sensitivity profile (adjoint method)
        sensi_profile = phi(:,1).*phi(:,i+1);

        %interpolate volume to sensitivity profile of shape (xdim, ydim, zdim)
        Sensitivity_profile(:,:,:,i) = griddata(cfg.node(:,1), cfg.node(:,2), cfg.node(:,3),sensi_profile, xi, yi, zi);

    
        %Normalize by the sum to get the density probability
        Sensitivity_profile(:,:,:,i) = Sensitivity_profile(:,:,:,i)/sum(Sensitivity_profile(:,:,:,i),"all");


        
        %Get fluence at fiber detector position
        
        
        %Select image around detector
        %Detector is at position [ceil(radius_fiber_det_mm)+1 , ceil(radius_fiber_det_mm)+1]
        Phi_around_det = Phi_volume(cfg.detpos(i,2)+1 - ceil(cfg.radius_fiber_det_mm) : cfg.detpos(i,2)+1 + ceil(cfg.radius_fiber_det_mm) , ...
                                    cfg.detpos(i,1)+1 - ceil(cfg.radius_fiber_det_mm) : cfg.detpos(i,1)+1 + ceil(cfg.radius_fiber_det_mm),1);
                                     
	    dPhi_dz_around_det = Fz(cfg.detpos(i,2)+1 - ceil(cfg.radius_fiber_det_mm) : cfg.detpos(i,2)+1 + ceil(cfg.radius_fiber_det_mm) , ...
                                cfg.detpos(i,1)+1 - ceil(cfg.radius_fiber_det_mm) : cfg.detpos(i,1)+1 + ceil(cfg.radius_fiber_det_mm),1);
                                                 
                                                 
 	    %Calculate diffuse reflectance
 	    DR = calculate_Diffuse_Reflectance_from_Fluence(Phi_around_det, dPhi_dz_around_det, 1, optical_prop.n_skin, optical_prop.mua_skin, (optical_prop.mu_s_skin)*(1-optical_prop.g_skin) );
        
        %Interpolate image around detector to match reso_detector_mm
        Rin  = imref2d(size(DR), [0 2*ceil(cfg.radius_fiber_det_mm)], [0 2*ceil(cfg.radius_fiber_det_mm)]); 
        Rout = imref2d([dim_interp_img dim_interp_img],   [0 2*ceil(cfg.radius_fiber_det_mm)], [0 2*ceil(cfg.radius_fiber_det_mm)]);
        
        tform = affine2d(eye(3)); % identity (no geometric change, only resampling)
        J = imwarp(DR, Rin, tform, ...
                   'OutputView', Rout, ...
                   'InterpolationMethod', 'cubic');
        
        % pixel center
        cx = (dim_interp_img+1)/2;  
        cy = (dim_interp_img+1)/2;
    
        %Mask on detector fiber position
        radius_fiber_det_px = cfg.radius_fiber_det_mm / Rout.PixelExtentInWorldX;
        [x, y] = meshgrid(1:dim_interp_img, 1:dim_interp_img);
        mask = (x - cx).^2 + (y - cy).^2 <= radius_fiber_det_px^2;
            
    
        DR_at_fiber_detector(end+1) = sum(J(mask));

    end
end
