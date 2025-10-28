function [data_around_det] = get_data_around_detector(data_surf, detpos, radius_fiber_det_mm)
%Get data around detector
% Data is a surface 2D array
    data_around_det = data_surf(detpos(2)+1 - ceil(radius_fiber_det_mm) : detpos(2)+1 + ceil(radius_fiber_det_mm) , ...
                                detpos(1)+1 - ceil(radius_fiber_det_mm) : detpos(1)+1 + ceil(radius_fiber_det_mm));
end