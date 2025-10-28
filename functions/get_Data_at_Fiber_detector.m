function [Data_at_fiber_detector] = get_Data_at_Fiber_detector(Data_around_detector, reso_detector_mm, radius_fiber_det_mm)
        % Interpolate image around detector to match reso_detector_mm and
        % get values summed at detector fiber position
        % Data_around_detector: Image of data around det (resolution 1mm)


        % Suppose A is an image of size Ny x Nx
        [Ny, Nx] = size(Data_around_detector);

        % Define original grid (resolution 1mm)
        [x, y] = meshgrid(1:Nx, 1:Ny);

        % Define new grid (resolution cfg.reso_detector_mm )
        [xq, yq] = meshgrid(linspace(1, Nx, ceil(Nx/reso_detector_mm)), ...
                            linspace(1, Ny, ceil(Ny/reso_detector_mm)));
        
        % Interpolate DR image to mach reslution cfg.reso_detector_m
        interp_Data = interp2(x, y, Data_around_detector, xq, yq, 'cubic');

        % Center of the fiber optics
        cx = (size(interp_Data,1)+1)/2;  
        cy = (size(interp_Data,2)+1)/2;

        %Mask on detector fiber position
        radius_fiber_det_px = ceil(radius_fiber_det_mm/reso_detector_mm);
        [x, y] = meshgrid(1:ceil(Nx/reso_detector_mm), ...
                          1:ceil(Ny/reso_detector_mm));

        mask = (x - cx).^2 + (y - cy).^2 <= radius_fiber_det_px^2;
        
        % Sum diffuse reflectance values at the level of the fiber detector
        Data_at_fiber_detector = sum(interp_Data(mask));
end