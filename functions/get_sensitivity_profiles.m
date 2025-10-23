function [Sensitivity_profile, phi_detval, Fluence_at_fiber_detector, DR_at_fiber_detector] = get_sensitivity_profiles(cfg, optical_prop)
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
    DR_at_fiber_detector = [];
    Fluence_at_fiber_detector = [];

    % Create volume of fluence
    [xi, yi, zi] = meshgrid(0.5:cfg.xdim_mm-0.5, 0.5:cfg.ydim_mm-0.5, 0.5:cfg.zdim_mm-0.5);
    Phi_volume = griddata(cfg.node(:,1), cfg.node(:,2), cfg.node(:,3), phi(:,1), xi, yi, zi);

    %Initi sensitiviy profile
    Sensitivity_profile = zeros(length(xi), length(yi), length(zi), length(cfg.detectors_SD_mm));

    %Calculate derivative of fluence along z axis (resolution is 1mm3)
    [Fx,Fy,Fz] = gradient(Phi_volume,1,1,1);
    


    %Calculate sensitivity profile per detector
    for i=1:length(cfg.detectors_SD_mm)

        %Sensitivity profile (adjoint method)
        sensi_profile = phi(:,1).*phi(:,i+1);

        %interpolate volume to sensitivity profile of shape (xdim, ydim, zdim)
        Sensitivity_profile(:,:,:,i) = griddata(cfg.node(:,1), cfg.node(:,2), cfg.node(:,3),sensi_profile, xi, yi, zi);

    
        %Normalize by the sum to get the density probability
        Sensitivity_profile(:,:,:,i) = Sensitivity_profile(:,:,:,i)/sum(Sensitivity_profile(:,:,:,i),"all");


        
        %Select image around detector (reso is 1mm)
        %Detector is at position [ceil(radius_fiber_det_mm)+1 , ceil(radius_fiber_det_mm)+1]
        Phi_around_det = get_data_around_detector(Phi_volume(:,:,1), cfg.detpos(i,:), cfg.radius_fiber_det_mm);
        dPhi_dz_around_det = get_data_around_detector(Fz(:,:,1), cfg.detpos(i,:), cfg.radius_fiber_det_mm);
                                                 
                                                 
 	    %Calculate diffuse reflectance
 	    DR = calculate_Diffuse_Reflectance_from_Fluence(Phi_around_det, dPhi_dz_around_det, 1, optical_prop.n_skin, optical_prop.mua_skin, (optical_prop.mu_s_skin)*(1-optical_prop.g_skin) );
        
        % Sum diffuse reflectance values at the level of the fiber detector
        DR_at_fiber_detector(end+1) = get_Data_at_Fiber_detector(DR, cfg.reso_detector_mm, cfg.radius_fiber_det_mm);
   
        % Sum fluence at detector value
        Fluence_at_fiber_detector(end+1) = get_Data_at_Fiber_detector(Phi_around_det, cfg.reso_detector_mm, cfg.radius_fiber_det_mm);
   

    end
end
