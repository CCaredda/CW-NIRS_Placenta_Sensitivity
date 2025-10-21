function R = fresnel_mean_reflectance(thetai, n1, n2)
	% FRESNEL_MEAN_REFLECTANCE computes the mean Fresnel reflection coefficient
	% for unpolarized light at a planar boundary.
	%
	% Inputs:
	%   thetai - incidence angle in radians (relative to normal)
	%   n1     - refractive index of incident medium
	%   n2     - refractive index of transmission medium
	%
	% Output:
	%   R      - mean reflection coefficient (0..1)
	%
	% Example:
	%   R = fresnel_mean_reflectance(pi/6, 1.0, 1.4);

	% Snell's law
	sin_theta_t = n1/n2 * sin(thetai);

	% Handle total internal reflection
	TIR = abs(sin_theta_t) > 1;
	sin_theta_t(TIR) = NaN; % will set R=1 later

	theta_t = asin(sin_theta_t);

	% Cosines
	cos_theta_i = cos(thetai);
	cos_theta_t = cos(theta_t);

	% Fresnel coefficients
	Rs = ((n1.*cos_theta_i - n2.*cos_theta_t)./(n1.*cos_theta_i + n2.*cos_theta_t)).^2;
	Rp = ((n1.*cos_theta_t - n2.*cos_theta_i)./(n1.*cos_theta_t + n2.*cos_theta_i)).^2;

	% Mean reflection coefficient for unpolarized light
	R = 0.5*(Rs + Rp);

	% Set reflection to 1 for total internal reflection
	R(TIR) = 1;

end

