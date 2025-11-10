function [DR_at_fiber_detector,Fluence_at_fiber_detector] = get_Diffuse_reflectance_from_Fluence_MCX(Phi_volume, detpos, radius_fiber_det_mm, reso_detector_mm, optical_prop)

    %Init output
    DR_at_fiber_detector = zeros(size(detpos,1),1);
    Fluence_at_fiber_detector = zeros(size(DR_at_fiber_detector));

    %Calculate derivative of fluence along z axis (resolution is 1mm3)
    [Fx,Fy,Fz] = gradient(Phi_volume(:,:,2:end),1,1,1);
    


    %Calculate sensitivity profile per detector
    for i=1:size(detpos,1)
       
        %Select image around detector (reso is 1mm)
        %Detector is at position [ceil(radius_fiber_det_mm)+1 , ceil(radius_fiber_det_mm)+1]
        % Keep second slice for Phi (inside skin) because first slice is
        % air
        Phi_around_det = get_data_around_detector(Phi_volume(:,:,2), detpos(i,:), radius_fiber_det_mm);
         % Keep first slice for dPhi_dz (inside skin) because first slice
         % has not been considered in gradient calculation
        dPhi_dz_around_det = get_data_around_detector(Fz(:,:,1), detpos(i,:), radius_fiber_det_mm);
                                                 
                                                 
 	    %Calculate diffuse reflectance
 	    DR = calculate_Diffuse_Reflectance_from_Fluence(Phi_around_det, dPhi_dz_around_det, 1, optical_prop.n_skin, optical_prop.mua_skin, (optical_prop.mu_s_skin)*(1-optical_prop.g_skin) );
        
        % Sum diffuse reflectance values at the level of the fiber detector
        DR_at_fiber_detector(i) = get_Data_at_Fiber_detector(DR, reso_detector_mm, radius_fiber_det_mm);
   
        % Sum fluence at detector value
        Fluence_at_fiber_detector(i) = get_Data_at_Fiber_detector(Phi_around_det, reso_detector_mm, radius_fiber_det_mm);
   

    end
end
