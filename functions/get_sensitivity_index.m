function [Sensitivity_indexes] = get_sensitivity_index(cfg, sensitivity_profile, thickness_layers_mm)
    
    %Inputs
    %cfg: redbird structure that containes mesh info
    % Sensitivity_profile: Sensitivity profile shape(nodes, SD separation)
    % thickness_layers_mm: thickness layers (skin, adipose tissue, muscle)

    tissue_idx = zeros(length(thickness_layers_mm)+1,2);
    tissue_idx(1,:) = [1 1+thickness_layers_mm(1)];
    for i=2:length(thickness_layers_mm)
        s = 1+thickness_layers_mm(i-1);
        tissue_idx(i,:) = [s s+thickness_layers_mm(i)];
    end
    tissue_idx(end,:) = [1+thickness_layers_mm(end) cfg.zdim_mm];


    Sensitivity_indexes = zeros(length(cfg.detectors_SD_mm),size(tissue_idx,1));

    for sd=1:length(cfg.detectors_SD_mm)

        %Get sensitivity indexes
        for m=1:size(tissue_idx,1)
            Sensitivity_indexes(sd,m) = sum(sensitivity_profile(:, :, tissue_idx(m,1):tissue_idx(m,2), sd),"all");
        end
    end

end