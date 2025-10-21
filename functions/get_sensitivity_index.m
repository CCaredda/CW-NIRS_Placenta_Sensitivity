function [Sensitivity_indexes] = get_sensitivity_index(cfg, sensitivity_profile, thickness_layers_mm)
    
    %Inputs
    %cfg: redbird structure that containes mesh info
    % Sensitivity_profile: Sensitivity profile shape(nodes, SD separation)
    % thickness_layers_mm: thickness layers (skin, adipose tissue, muscle)

    slice_tissue = 4*ones(cfg.zdim_mm,1);
    slice_tissue(1:thickness_layers_mm(1)) = 1;
    slice_tissue(thickness_layers_mm(1)+1:thickness_layers_mm(1)+thickness_layers_mm(2)) = 2;
    slice_tissue(thickness_layers_mm(1)+thickness_layers_mm(2)+1:thickness_layers_mm(1)+thickness_layers_mm(2)+thickness_layers_mm(3)) = 3;



    Sensitivity_indexes = zeros(length(cfg.detectors_SD_mm),4);

    for sd=1:length(cfg.detectors_SD_mm)

        %Get sensitivity indexes
        for m=1:4
            % Find indices of tissue m
            indices = find(slice_tissue == m);
            % Compute sensitivity
            Sensitivity_indexes(sd,m) = sum(sensitivity_profile(:, :, indices,sd),"all");
        end
    end

end