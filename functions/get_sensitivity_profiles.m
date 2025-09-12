function [Sensitivity_profile, Diffuse_reflectance] = get_sensitivity_profiles(cfg, optical_prop)
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

    %Calculate sensitivity profile per detector, size(nodes,detectors)
    Sensitivity_profile = zeros(size(phi,1),length(cfg.detectors_SD_mm));
    for i=1:length(cfg.detectors_SD_mm)
        %Sensitivity profile (adjoint method)
        Sensitivity_profile(:,i) = phi(:,1).*phi(:,i+1);
    
        %Normalize by the sum to get the density probability
        Sensitivity_profile(:,i) = Sensitivity_profile(:,i)/sum(Sensitivity_profile(:,i),"all");
    end
end