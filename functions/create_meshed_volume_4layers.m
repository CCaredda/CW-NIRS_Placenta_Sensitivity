function [cfg] = create_meshed_volume_4layers(zdim_mm, xdim_mm, ydim_mm, thickness_layers_mm, max_vol_mesh, detectors_SD_mm, display)
    %Create layers mesh volume
    %Input:
    % volume_depth_mm, volume_width_mm, volume_length_mm: dimension of the volume in mm
    % thickness_layers_mm: thickness of the 3 first layer in mm (the thickness of
    % the fourth layer is calculated with volume_depth_mm
    % max_vol_mesh: Maximum tetrahedral element volume (in cubic voxel-length unit).
    % For example, setting 𝑉max to 4 means no tetrahedral element in the
    % generated mesh can exceed 4 cubic voxels in volume. For example  max_vol_mesh= [0.1; 0.1; 1; 1000];
    % detectors_SD_mm: source detector separation in mm
    % display: 0 or 1 fr displaying the volume
    
    %Source detector separation
    cfg.detectors_SD_mm  = detectors_SD_mm;

    %Thickness of the layers
    N_layers = 4;
    thickness_slices_mm = [0];
    for i=1:N_layers-1
        thickness_slices_mm(end+1) = sum(thickness_layers_mm(1:i));
    end
    thickness_slices_mm(end+1) = zdim_mm;

    
    %Define volume
    [cfg.node,cfg.face,c0]=latticegrid([0 xdim_mm],...
                                       [0 ydim_mm],...
                                       thickness_slices_mm); % Create the four-layered mesh surface;
    fc2=cell2mat(cfg.face);
    cfg.face=[fc2(:,[1 2 3]); fc2(:,[1 3 4])];
    
    % per-layer target tetra volumes (example numbers)
    c0(:,4) = max_vol_mesh;
    
    %Convert surf to node
    [cfg.node,cfg.elem]=surf2mesh(cfg.node,cfg.face,[],[],1,[],c0);
    
    %Get label per node
    labels = cfg.elem(:,5);
    
    
    % % Get unique element idx and corresponding labels that are related to the nodes
    % cfg.unique_elems = [];
    % cfg.labels_per_node = [];
    % for i=1:size(cfg.elem,2)
    %     [val, idx, ~] = unique(cfg.elem(:,i));
    %     cfg.unique_elems = [cfg.unique_elems;val];
    %     cfg.labels_per_node = [cfg.labels_per_node;labels(idx)];
    % end
    % [cfg.unique_elems,unique_elem_idx,~] = unique(cfg.unique_elems);
    % cfg.labels_per_node = cfg.labels_per_node(unique_elem_idx);
    
    
    % Display model
    if display
        figure;  plottetra(cfg.node,[cfg.elem labels]);
    end


    %Store volume dim
    cfg.xdim_mm = xdim_mm;
    cfg.ydim_mm = ydim_mm;
    cfg.zdim_mm = zdim_mm;
    



    %Info detectors
    %detectors_SD_mm = [30, 40, 50];
    
    %Light source
    src_pos = [(xdim_mm/2-1) (ydim_mm - detectors_SD_mm(end))/2-1 0];
    cfg.srcdir=[0 0 1];
    cfg.srcpos = src_pos;
    cfg.srctype = 'pencil';
    
    %Detectors
    %Do not model the fiber optics (it is scaled with the calibration procedure)
    cfg.detdir = [0 0 1];
    det_pos = zeros(length(detectors_SD_mm),3);
    for i=1:length(detectors_SD_mm)
        det_pos(i,:) = [src_pos(1) src_pos(2)+detectors_SD_mm(i) 0]; % No radius!!
    end
    cfg.detpos = det_pos;
    
    
    % modulation angular frequency
    cfg.omega = 0;

end