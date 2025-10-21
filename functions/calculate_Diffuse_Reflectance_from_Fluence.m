function [DR] = calculate_Diffuse_Reflectance_from_Fluence(Phi_surf,dPhi_dz_surf,n_air,n_skin,mua, musp )
    %Phi_surf: Image of Fluence value at z=0 
    % dPhi_dz_surf: Image of derivative along z axis of Fluence value at z=0 
    
    
    d_theta = 0.05;
    theta = 0:d_theta:2*pi;
    R = fresnel_mean_reflectance(theta,n_skin,n_air);
    
    

    D = 1/(3*(mua+musp));

    DR = zeros(size(Phi_surf));
    
    for x = 1:size(Phi_surf,1)
        for y = 1:size(Phi_surf,2)
            DR(x,y) = sum(d_theta*(1-R)*(1/(4*pi)) .* (Phi_surf(x,y) + 3*D*dPhi_dz_surf(x,y)*cos(theta)) .* cos(theta));
        end
    end
end
