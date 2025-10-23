function [DR] = calculate_Diffuse_Reflectance_from_Fluence(F,dF_dz, n_air, n_skin, mua, musp )
    %F: Image of Fluence value at z=0 
    % dF_dz: Image of derivative along z axis of Fluence value at z=0 
    %Eq 7 Alwin Kienle and Michael S. Patterson, "Improved solutions of the steady-state and the time-resolved diffusion equations for reflectance from a semi-infinite turbid medium," J. Opt. Soc. Am. A 14, 246-254 (1997)
    
    d_theta = 0.05;
    theta = 0:d_theta:pi/2;
   
    %Calculate Fresnel reflection coefficient 
    R = fresnel_mean_reflectance(theta,n_skin,n_air);
    
    

    D = 1/(3*(mua+musp));

    DR = zeros(size(F));
    
    %Calculate diffuse reflectance
    for x = 1:size(F,1)
        for y = 1:size(F,2)
            DR(x,y) = sum((1-R)*(1/(4*pi)) .* (F(x,y) + 3*D*dF_dz(x,y)*cos(theta)) .* cos(theta) .* sin(theta) *d_theta)*2*pi;
        end
    end
end
