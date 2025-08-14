clear
close all
clc

cluster = 0;

% Add path
if cluster ==1
    addpath('~/private/redbird/matlab');
    addpath('~/private/iso2mesh');
else
    addpath('~/Soft/redbird/matlab');
    addpath('~/Soft/iso2mesh');
end

addpath('../functions');



%Wavelength
Lambdas = 780;
% Lambdas = 840;
% Lambdas = 890;
% Lambdas = [780, 840, 890 ];


%Arrays thickness
skin_thickness_array_mm = [1,2,3,4,5];
adipose_thickness_array_mm = [1,2,4,8,16];
muscle_thickness_array_mm = [3,7,10,13,20,25];

%Saturation array
SatO2_muscle_array = 0.4:0.2:0.8;
SatO2_placenta_array = 0.4:0.2:0.8;


%HbT array
C_HbT_muscle_array = [15,25,35,50]*1e-6;
C_HbT_placenta_array = [15,25,35,50]*1e-6;


%Volume fraction of melanosome according to the color tones
% Modeling and Verification of Melanin Concentration on Human Skin Type
f_melanosome = [0.0255 0.155 0.305];


%Info volume
volume_width_mm = 70;
volume_length_mm = 110;
volume_depth_mm = 40;

%Info detectors
detectors_SD_mm = [30, 40, 50];

%Light source
src_pos = [(volume_width_mm/2-1) (volume_length_mm - detectors_SD_mm(end))/2-1 0];
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

%Save constant
if ~exist('output_lookup_table', 'dir')
    mkdir('output_lookup_table')
end
save('output_lookup_table/cst.mat','volume_width_mm','volume_length_mm','volume_depth_mm','src_pos','det_pos','Lambdas','SatO2_placenta_array','SatO2_muscle_array','skin_thickness_array_mm','adipose_thickness_array_mm','muscle_thickness_array_mm','f_melanosome');









%Define Volume
thickness_skin_in_mm = skin_thickness_array_mm(1);
thickness_adipose_in_mm = adipose_thickness_array_mm(1);
thickness_muscle_in_mm = muscle_thickness_array_mm(1);

%Define volume
tic
[cfg.node,cfg.face,c0]=latticegrid([0 volume_width_mm],...
                                   [0 volume_length_mm],...
                                   [0 thickness_skin_in_mm thickness_skin_in_mm+thickness_adipose_in_mm thickness_skin_in_mm+thickness_adipose_in_mm+thickness_muscle_in_mm volume_depth_mm]); % Create the four-layered mesh surface;
fc2=cell2mat(cfg.face);
cfg.face=[fc2(:,[1 2 3]); fc2(:,[1 3 4])];

% per-layer target tetra volumes (example numbers)
c0(:,4) = [0.1; 0.1; 1; 1000];

%Convert surf to node
[cfg.node,cfg.elem]=surf2mesh(cfg.node,cfg.face,[],[],1,[],c0);
%Get label per node
labels = cfg.elem(:,5);


% Get unique element idx and corresponding labels that are related to the nodes
unique_elems = [];
labels_per_node = [];
for i=1:size(cfg.elem,2)
    [val, idx, ~] = unique(cfg.elem(:,i));
    unique_elems = [unique_elems;val];
    labels_per_node = [labels_per_node;labels(idx)];
end
[unique_elems,unique_elem_idx,~] = unique(unique_elems);
labels_per_node = labels_per_node(unique_elem_idx);

toc
%TEMP Display model
% figure;  plottetra(cfg.node,[cfg.elem labels]);




%% Test with one wavelength
% 2.419216 * 432 * 150 /(60*60*24) (s) -> 1.8 days


tic
% disp('Define optical properties started');
%Temp
C_HbT_muscle = C_HbT_muscle_array(1);
C_HbT_placenta = C_HbT_placenta_array(1);
SatO2_muscle = SatO2_muscle_array(1);
SatO2_placenta = SatO2_placenta_array(1);
f_mel = f_melanosome(1);

%Define optical properties for all wavelengths
cfg.prop = containers.Map();
% Calculate optical properties
optical_prop = process_optical_properties_skin_Fat_muscle_placenta(Lambdas,f_mel,SatO2_muscle, SatO2_placenta,C_HbT_muscle,C_HbT_placenta);

% Set optical properties % [mua,mus,g,n]
% 0: Air
% 1: Skin
% 2: Adipose tissue
% 3: Muscle
% 4: Placenta
cfg.prop=[0 0 1 1; ...
optical_prop.mua_skin optical_prop.mu_s_skin optical_prop.g_skin optical_prop.n_skin ; ...
optical_prop.mua_adipose optical_prop.mu_s_adipose optical_prop.g_adipose optical_prop.n_adipose ; ...
optical_prop.mua_muscle optical_prop.mu_s_muscle optical_prop.g_muscle optical_prop.n_muscle ; ...
optical_prop.mua_placenta optical_prop.mu_s_placenta optical_prop.g_placenta optical_prop.n_placenta];








% %Loop on HbT muscle
% for C_HbT_muscle = C_HbT_muscle_array
%     for C_HbT_placenta = C_HbT_placenta_array
%         % Loop on SatO2 muscle
%         for SatO2_muscle = SatO2_muscle_array
%             % Loop on SatO2 placenta
%             for SatO2_placenta = SatO2_placenta_array
%                 %Loop on color tones
%                 for f_mel = f_melanosome
%                     %Define optical properties for all wavelengths
%                     cfg.prop = containers.Map();
%                     % Calculate optical properties
%                     optical_prop = process_optical_properties_skin_Fat_muscle_placenta(Lambdas,f_mel,SatO2_muscle, SatO2_placenta,C_HbT_muscle,C_HbT_placenta);
% 
%                     % Set optical properties % [mua,mus,g,n]
%                     % 0: Air
%                     % 1: Skin
%                     % 2: Adipose tissue
%                     % 3: Muscle
%                     % 4: Placenta
%                     cfg.prop=[0 0 1 1; ...
%                     optical_prop.mua_skin optical_prop.mu_s_skin optical_prop.g_skin optical_prop.n_skin ; ...
%                     optical_prop.mua_adipose optical_prop.mu_s_adipose optical_prop.g_adipose optical_prop.n_adipose ; ...
%                     optical_prop.mua_muscle optical_prop.mu_s_muscle optical_prop.g_muscle optical_prop.n_muscle ; ...
%                     optical_prop.mua_placenta optical_prop.mu_s_placenta optical_prop.g_placenta optical_prop.n_placenta];
% 
% 
%                 end
%             end
%         end
%     end
% end



%Prepare all missing values in cfg struct
% disp('Forward solver...')
cfg = rbmeshprep(cfg);


%Forward solver
[Diffuse_reflectance, phi] = rbrun(cfg);
%Set negative values to 0
phi(phi<0) = 0;



% phi(phi<0)=0;
%phi: the full volumetric forward solution computed at all wavelengths at every source and detector
%Temp
% close all
% figure, plotmesh([cfg.node, log10(abs(phi(:,1)))], cfg.elem, strcat('y=',num2str(volume_length_mm/2)), 'facecolor', 'interp', 'linestyle', 'none');
% figure, plotmesh([cfg.node, log10(abs(phi(:,4)))], cfg.elem, strcat('y=',num2str(volume_length_mm/2)), 'facecolor', 'interp', 'linestyle', 'none');
% figure, plotmesh([cfg.node, log10(abs(phi(:,1)))], cfg.elem, strcat('x=',num2str(volume_width_mm/2)), 'facecolor', 'interp', 'linestyle', 'none');
% figure, plotmesh([cfg.node, log10(abs(phi(:,4)))], cfg.elem, strcat('x=',num2str(volume_width_mm/2)), 'facecolor', 'interp', 'linestyle', 'none');
figure, plotmesh([cfg.node, log10(abs(phi(:,1).*phi(:,4)))], cfg.elem, strcat('x=',num2str(volume_width_mm/2)), 'facecolor', 'interp', 'linestyle', 'none');



    

%Calculate sensitivity profile per detector, size(nodes,detectors)
Sensitivity_profile = zeros(size(phi,1),length(detectors_SD_mm));
for i=1:length(detectors_SD_mm)
    %Sensitivity profile (adjoint method)
    Sensitivity_profile(:,i) = phi(:,1).*phi(:,i+1);

    %Normalize by the sum to get the density probability
    Sensitivity_profile(:,i) = Sensitivity_profile(:,i)/sum(Sensitivity_profile(:,i),"all");
end


%Get sensitivity indexes
idx_labels = unique(labels_per_node);
Sensitivity_indexes = zeros(length(detectors_SD_mm),length(idx_labels));
for m=1:length(idx_labels)
    % Find indices of unique_elems related to tissue m
    idx_elem = find(labels_per_node == idx_labels(m));

    %Get elem of the tissue.
    %It corresponded of the index of the tetrahedron in cfg.node
    elem_tissue = unique_elems(idx_elem,:);

    % Compute sensitivity
    Sensitivity_indexes(:,m) = sum(Sensitivity_profile(elem_tissue,:),1);
end

%Save outputs
output_name = strcat('output_lookup_table/out_St_muscle_',num2str(SatO2_muscle),'_St_placenta_',num2str(SatO2_placenta),'_Thick_skin_',num2str(thickness_skin_in_mm),'_Thick_adipose_',num2str(thickness_adipose_in_mm),'_Thick_muscle_',num2str(thickness_muscle_in_mm),'f_mel',num2str(f_mel),'_HbT_muscle_umol_',num2str(C_HbT_muscle*1e6),'_HbT_placenta_umol_',num2str(C_HbT_placenta*1e6),'.mat');
save(output_name,'Diffuse_reflectance','Sensitivity_indexes');

toc





























%% Test with multispectral (6 or 2 wavelenths)

% 11.409799 * 432 * 150 /(60*60*24) (s) -> 8.5 days (6 wavelength)
% 7.230005 * 432 * 150 /(60*60*24) -> 5.4 days (2 wavelength)
   

% tic
% % disp('Define optical properties started');
% %Temp
% C_HbT_muscle = C_HbT_muscle_array(1);
% C_HbT_placenta = C_HbT_placenta_array(1);
% SatO2_muscle = SatO2_muscle_array(1);
% SatO2_placenta = SatO2_placenta_array(1);
% f_mel = f_melanosome(1);
% 
% %Define optical properties for all wavelengths
% cfg.prop = containers.Map();
% % Calculate optical properties
% optical_prop = process_optical_properties_skin_Fat_muscle_placenta(Lambdas,f_mel,SatO2_muscle, SatO2_placenta,C_HbT_muscle,C_HbT_placenta);
% 
% % Set optical properties % [mua,mus,g,n]
% % 0: Air
% % 1: Skin
% % 2: Adipose tissue
% % 3: Muscle
% % 4: Placenta
% for l=1:length(Lambdas)
%     cfg.prop(num2str(Lambdas(l)))=[0 0 1 1; ...
%     optical_prop.mua_skin(l) optical_prop.mu_s_skin(l) optical_prop.g_skin optical_prop.n_skin ; ...
%     optical_prop.mua_adipose(l) optical_prop.mu_s_adipose(l) optical_prop.g_adipose optical_prop.n_adipose ; ...
%     optical_prop.mua_muscle(l) optical_prop.mu_s_muscle(l) optical_prop.g_muscle optical_prop.n_muscle ; ...
%     optical_prop.mua_placenta(l) optical_prop.mu_s_placenta(l) optical_prop.g_placenta optical_prop.n_placenta];
% end
% 
% 
% 
% 
% 
% 
% 
% % %Loop on HbT muscle
% % for C_HbT_muscle = C_HbT_muscle_array
% %     for C_HbT_placenta = C_HbT_placenta_array
% %         % Loop on SatO2 muscle
% %         for SatO2_muscle = SatO2_muscle_array
% %             % Loop on SatO2 placenta
% %             for SatO2_placenta = SatO2_placenta_array
% %                 %Loop on color tones
% %                 for f_mel = f_melanosome
% %                     %Define optical properties for all wavelengths
% %                     cfg.prop = containers.Map();
% %                     % Calculate optical properties
% %                     optical_prop = process_optical_properties_skin_Fat_muscle_placenta(Lambdas,f_mel,SatO2_muscle, SatO2_placenta,C_HbT_muscle,C_HbT_placenta);
% % 
% %                     % Set optical properties % [mua,mus,g,n]
% %                     % 0: Air
% %                     % 1: Skin
% %                     % 2: Adipose tissue
% %                     % 3: Muscle
% %                     % 4: Placenta
% %                     for l=1:length(Lambdas)
% %                         cfg.prop(num2str(Lambdas(l)))=[0 0 1 1; ...
% %                         optical_prop.mua_skin(l) optical_prop.mu_s_skin(l) optical_prop.g_skin optical_prop.n_skin ; ...
% %                         optical_prop.mua_adipose(l) optical_prop.mu_s_adipose(l) optical_prop.g_adipose optical_prop.n_adipose ; ...
% %                         optical_prop.mua_muscle(l) optical_prop.mu_s_muscle(l) optical_prop.g_muscle optical_prop.n_muscle ; ...
% %                         optical_prop.mua_placenta(l) optical_prop.mu_s_placenta(l) optical_prop.g_placenta optical_prop.n_placenta];
% %                     end
% % 
% %                 end
% %             end
% %         end
% %     end
% % end
% 
% 
% 
% %Prepare all missing values in cfg struct
% % disp('Forward solver...')
% cfg = rbmeshprep(cfg);
% 
% 
% %Forward solver
% [dr_w, phi_w] = rbrun(cfg);
% 
% 
% 
% 
% % phi(phi<0)=0;
% %phi: the full volumetric forward solution computed at all wavelengths at every source and detector
% %Temp
% % val = phi_w(num2str(Lambdas(1)));
% % close all
% % figure, plotmesh([cfg.node, log10(abs(val(:,1)))], cfg.elem, strcat('y=',num2str(volume_length_mm/2)), 'facecolor', 'interp', 'linestyle', 'none');
% % figure, plotmesh([cfg.node, log10(abs(val(:,4)))], cfg.elem, strcat('y=',num2str(volume_length_mm/2)), 'facecolor', 'interp', 'linestyle', 'none');
% % figure, plotmesh([cfg.node, log10(abs(val(:,1)))], cfg.elem, strcat('x=',num2str(volume_width_mm/2)), 'facecolor', 'interp', 'linestyle', 'none');
% % figure, plotmesh([cfg.node, log10(abs(val(:,4)))], cfg.elem, strcat('x=',num2str(volume_width_mm/2)), 'facecolor', 'interp', 'linestyle', 'none');
% %figure, plotmesh([cfg.node, log10(abs(val(:,1).*val(:,4)))], cfg.elem, strcat('x=',num2str(volume_width_mm/2)), 'facecolor', 'interp', 'linestyle', 'none');
% 
% 
% for l=Lambdas
%     %Get phi and diffuse reflectance values at specific wavelength
%     phi = phi_w(num2str(l));
%     Diffuse_reflectance = dr_w(num2str(l));
% 
%     %Set negative values to 0
%     phi(phi<0) = 0;
% 
%     %Calculate sensitivity profile per detector, size(nodes,detectors)
%     Sensitivity_profile = zeros(size(phi,1),length(detectors_SD_mm));
%     for i=1:length(detectors_SD_mm)
%         %Sensitivity profile (adjoint method)
%         Sensitivity_profile(:,i) = phi(:,1).*phi(:,i+1);
% 
%         %Normalize by the sum to get the density probability
%         Sensitivity_profile(:,i) = Sensitivity_profile(:,i)/sum(Sensitivity_profile(:,i),"all");
%     end
% 
% 
%     %Get sensitivity indexes
%     idx_labels = unique(labels_per_node);
%     Sensitivity_indexes = zeros(length(detectors_SD_mm),length(idx_labels));
%     for m=1:length(idx_labels)
%         % Find indices of unique_elems related to tissue m
%         idx_elem = find(labels_per_node == idx_labels(m));
% 
%         %Get elem of the tissue.
%         %It corresponded of the index of the tetrahedron in cfg.node
%         elem_tissue = unique_elems(idx_elem,:);
% 
%         % Compute sensitivity
%         Sensitivity_indexes(:,m) = sum(Sensitivity_profile(elem_tissue,:),1);
%     end
% 
%     %Save outputs
%     output_name = strcat('output_lookup_table/out_St_muscle_',num2str(SatO2_muscle),'_St_placenta_',num2str(SatO2_placenta),'_Thick_skin_',num2str(thickness_skin_in_mm),'_Thick_adipose_',num2str(thickness_adipose_in_mm),'_Thick_muscle_',num2str(thickness_muscle_in_mm),'f_mel',num2str(f_mel),'_HbT_muscle_umol_',num2str(C_HbT_muscle*1e6),'_HbT_placenta_umol_',num2str(C_HbT_placenta*1e6),'_wavelength_',num2str(l),'.mat');
%     save(output_name,'Diffuse_reflectance','Sensitivity_indexes');
% 
% end
% toc



                 
                    
                 
