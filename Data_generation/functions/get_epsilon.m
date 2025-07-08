function [eps_hb,eps_hb02,eps_oxCCO,eps_redCCO] = get_epsilon(lambda)


%Find lambda values into w
w = readmatrix('spectra/lambda.txt');
[sharedVals,idx_w] = intersect(w,lambda,'stable');

%Hb
eps = readmatrix('spectra/eps_Hb.txt');
eps_hb = eps(idx_w);
%HbO2
eps = readmatrix('spectra/eps_HbO2.txt');
eps_hb02 = eps(idx_w);
%oxCCO
eps = readmatrix('spectra/eps_oxCCO.txt');
eps_oxCCO = eps(idx_w);
%redCCO
eps = readmatrix('spectra/eps_redCCO.txt');
eps_redCCO = eps(idx_w);

% %% Read the molar absorption coefficients [M^-1*cm^-1] of HbO2 and HHb from external file:
% % hemo    = load('hemoglobine.mat'); %http://omlc.org/spectra/hemoglobin/summary.html
% filename='Molar extinction coeff HbO2 & HHb.txt'; % Name of the text file containing the molar extinction coefficients of oxidized CCO;
% data=dlmread(filename,'\t',2,0); % Read the file containing the molar extinction coefficients of oxidized CCO into a matrix;
% eps_Hb_HbO2_lambda = data(:,1);
% eps_HbO2 = data(:,2);
% eps_Hb = data(:,3);
% 
% 
% % eps_hb = interp1(eps_Hb_HbO2_lambda,eps_Hb,lambda,'linear');
% [xData, yData] = prepareCurveData( eps_Hb_HbO2_lambda,eps_Hb);
% ft = 'linearinterp';
% [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
% eps_hb=feval(fitresult,lambda);%cm^-1/(mol.l^-1)
% 
% % eps_hb02 = interp1(eps_Hb_HbO2_lambda,eps_HbO2,lambda,'linear');
% [xData, yData] = prepareCurveData(eps_Hb_HbO2_lambda,eps_HbO2);
% ft = 'linearinterp';
% [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
% eps_hb02=feval(fitresult,lambda);%cm^-1/(mol.l^-1)
% 
% 
% 
% 
% %% Read the molar extinction coefficients [mM^-1*cm^-1] of oxidized and reduced CCO from external file:
% filename='Molar extinction coeff oxCCO.txt'; % Name of the text file containing the molar extinction coefficients of oxidized CCO;
% eps_oxCCO=dlmread(filename,'\t',2,0); % Read the file containing the molar extinction coefficients of oxidized CCO into a matrix;
% eps_oxCCO(:,2)=eps_oxCCO(:,2).*1000; % Unit conversion of the molar absorption coefficients of oxidized CCO, from [mM^-1*cm^-1] in [M^-1*cm^-1];
% 
% filename='Molar extinction coeff redCCO.txt'; % Name of the text file containing the molar extinction coefficients of reduced CCO;
% eps_redCCO=dlmread(filename,'\t',2,0); % Read the file containing the molar extinction coefficients of reduced CCO into a matrix;
% eps_redCCO(:,2)=eps_redCCO(:,2).*1000; % Unit conversion of the molar absorption coefficients of reduced CCO, from [mM^-1*cm^-1] in [M^-1*cm^-1];
% 
% 
% %interpolate mua values to the given Lambda range
% % eps_oxCCO = interp1(eps_oxCCO(:,1),eps_oxCCO(:,2),lambda,'linear');
% [xData, yData] = prepareCurveData(eps_oxCCO(:,1),eps_oxCCO(:,2));
% ft = 'linearinterp';
% [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
% eps_oxCCO=feval(fitresult,lambda);%cm^-1/(mol.l^-1)
% 
% % eps_redCCO = interp1(eps_redCCO(:,1),eps_redCCO(:,2),lambda,'linear');
% [xData, yData] = prepareCurveData(eps_redCCO(:,1),eps_redCCO(:,2));
% ft = 'linearinterp';
% [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
% eps_redCCO=feval(fitresult,lambda);%cm^-1/(mol.l^-1)





end