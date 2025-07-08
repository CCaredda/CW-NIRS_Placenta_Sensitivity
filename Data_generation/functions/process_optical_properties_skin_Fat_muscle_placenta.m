function [opt] = process_optical_properties_skin_Fat_muscle_placenta(lambda,f_mel,SatO2_muscle,SatO2_placenta,C_HbT_muscle,C_HbT_placenta)
    %Compute optical properties for the different classes
    %INPUT:
    %Lambda: wavelength (nm)
    %Output:
    %optical_prop: vector of optical properties size(length(Lambdas),7,4))


    %Find lambda values into w
    w = readmatrix('spectra/lambda.txt');
    [sharedVals,idx_w] = intersect(w,lambda,'stable');
    
    mua = readmatrix('spectra/mua_H2O.txt');
    mu_a_H2O = mua(idx_w);
    
    
    % Read the absorption coefficients [cm^-1] of fat from external file:
    mua = readmatrix('spectra/mua_Fat.txt');
    mu_a_fat = mua(idx_w);

    % Get extinction coefficient (in cm-.mol-1.L)
    [eps_hb,eps_hb02,eps_oxCCO,eps_redCCO] = get_epsilon(lambda);



    %% Skin

    % Ratio epidermis and dermis: Measurement of Epidermis, Dermis, and Total Skin Thicknesses from Six Different Body Regions with a New Ethical Histometric Technique
    T_epi = 127.2; %µm
    T_dermi = 4550; %µm
    F_dermis = T_dermi/(T_dermi+T_epi);
    F_epi = T_epi/(T_dermi+T_epi);

    % Mean values for dermin, epidermis and strum corneus
    opt.g_skin=0.83; % Modeling and Verification of Melanin Concentration on Human Skin Type
    
    % Mean values for dermin, epidermis and strum corneus
    opt.n_skin=1.47; % Modeling and Verification of Melanin Concentration on Human Skin Type

    %Absorption coeff of melanosome (Jacques)
    mua_mel = 519*(lambda/500).^(-3); % in cm-1

    %Skin baseline (Iyad S. Saidi, 1992 Transcutaneous optical measurement of hyperbilirubinemia in neonates. Ph.D. dissertation, Rice University, Houston, TX, USA )
    mua_skinbaseline = (7.84*1e8)*(lambda.^(-3.255));  %cm-1
    
    %Absorption epidermis coefficient Jacques
    mua_epi = f_mel*mua_mel  +  (1 - f_mel)*mua_skinbaseline;

    %Absorption dermis coefficient Jacques
    F_B = 0.002; % Volume fraction of blood in dermis  JACQUES
    C_HbT = 4.7 *1e-6; % CHbT in dermis (Jacques)
    SatO_2 = 0.39; %jacques

    mua_derm = (1-F_B)*mua_skinbaseline + log(10)*C_HbT*(SatO_2*eps_hb02 + (1-SatO_2)*eps_hb);

    mua_skin = F_dermis*mua_derm + F_epi*mua_epi;
    opt.mua_skin = mua_skin*0.1; % convert into mm-1
    
    % Scattering coeffcient
    mu_s_reduced=46*((lambda./500).^(-1.421));  % Reduced scattering coefficient of the tissue [cm^-1];
    mu_s_skin=mu_s_reduced./(1-opt.g_skin); % Scattering coefficient of the tissue [cm^-1];
    opt.mu_s_skin = 0.1*mu_s_skin; % convert into mm-1


    %% Adipose tissue
    
    opt.g_adipose=0.9;
    opt.n_adipose=1.37;

    Water_adipose = 0.8;
    Fat_adipose = 0.2;
    C_HbO2_adipose = 0;
    C_Hb_adipose = 0;

    % Absorption coefficient
    mua_adipose = Water_adipose*mu_a_H2O + ...
    Fat_adipose*mu_a_fat + ...
    log(10)*C_Hb_adipose*eps_hb + ...
    log(10)*C_HbO2_adipose*eps_hb02; % Absorption coefficient of the tissue [cm^-1];

    opt.mua_adipose = 0.1*mua_adipose; % convert into mm-1

    % Scattering coeffcient
    mu_s_reduced=18.4*((lambda./500).^(-0.672));  % Reduced scattering coefficient of the tissue [cm^-1];
    mu_s_adipose=mu_s_reduced./(1-opt.g_adipose); % Scattering coefficient of the tissue [cm^-1];
    opt.mu_s_adipose = 0.1*mu_s_adipose; % convert into mm-1

    %% Muscle tissue

    opt.g_muscle=0.9; %  Jacques p.8 (Samhatan 2012)
    opt.n_muscle=1.37;

    
    Water_muscle = 0.76;
    Fat_muscle =0;
    C_HbT = C_HbT_muscle;
    C_HbO2_muscle = SatO2_muscle*C_HbT;
    C_Hb_muscle = (1-SatO2_muscle)*C_HbT;

    mua_muscle = Water_muscle*mu_a_H2O + ...
    Fat_muscle*mu_a_fat + ...
    log(10)*C_Hb_muscle*eps_hb + ...
    log(10)*C_HbO2_muscle*eps_hb02; % Absorption coefficient of the tissue [cm^-1];
    
    % Absorption coefficient
    opt.mua_muscle = 0.1*mua_muscle; % convert into mm-1

    % Recovering fetal signals transabdominally through interferometric near-infrared spectroscopy (iNIRS)
    %mua = 0.01 mm-1 à 850 nm

    % Scattering coefficient
    mu_s_reduced=13*((lambda./500).^(-0.926)); % Reduced scattering coefficient of the tissue [cm^-1];
    % mu_s_reduced=9.8*((lambda./500).^(-2.82)); % Reduced scattering coefficient of the tissue [cm^-1];
    mu_s_muscle=mu_s_reduced./(1-opt.g_muscle); % Scattering coefficient of the tissue [cm^-1];
    opt.mu_s_muscle = mu_s_muscle*0.1; % Convert in mm-1

    %% Placenta tissue

    opt.g_placenta=0.9; % Recovering fetal signals transabdominally through interferometric near-infrared spectroscopy (iNIRS)
    opt.n_placenta=1.4;  % Recovering fetal signals transabdominally through interferometric near-infrared spectroscopy (iNIRS)

    Water_placenta = 0.85;
    Fat_placenta =0;
    C_HbT = C_HbT_placenta; % Non-invasive monitoring of blood oxygenation in human placentas via concurrent diffuse optical spectroscopy and ultrasound imaging
    C_HbO2_placenta = SatO2_placenta*C_HbT;
    C_Hb_placenta = (1-SatO2_placenta)*C_HbT;

    mua_placenta = Water_placenta*mu_a_H2O + ...
    Fat_placenta*mu_a_fat + ...
    log(10)*C_Hb_placenta*eps_hb + ...
    log(10)*C_HbO2_placenta*eps_hb02; % Absorption coefficient of the tissue [cm^-1];

    % Absorption coefficient
    opt.mua_placenta = 0.1*mua_placenta; % convert into mm-1
    
    % Scattering coefficient
    mu_s_reduced = 1.6619*((lambda./500).^(-1.426)); % Reduced scattering coefficient of the tissue [mm^-1];
    opt.mu_s_placenta=mu_s_reduced./(1-opt.g_placenta); % Scattering coefficient of the tissue [mm^-1];

    

end
