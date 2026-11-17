function result = fit_corona(out,figure_out)
% Verify the corona fit, unit conversion and the selected model inputs.
% Open this file and click Run. Only base MATLAB is required.
% Optional: result=fit_corona('absolute/output/folder').
%
% Source: Zheng, Zhang and He (2015), Phys. Plasmas 22, 023501.
% DOI 10.1063/1.4907234; negative corona in Fig. 3, PDF page 3
% (printed 023501-2). Fig. 5, PDF page 4, supplies a cross-check.
% Coordinates below were digitized from PDF vector marker centres.
% They are not raw experimental data. This script does not parse the PDF.
% The paper's current unit is microampere/cm (line current), not A or A/m^2.
%
% Conductivity reference: Foo et al. (2012), JAP 111, 034102,
% DOI 10.1063/1.3680878, Sec. V: VHB 4905 kappa0=3.23e-14 S/m,
% E_B=40 MV/m. Matching an input to this reference is not specimen validation.

root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin<1
 out=fullfile(root,'corona_validation','results');
 if nargin<2,figure_out=fullfile(root,'corona_validation','figure');end
elseif nargin<2,figure_out=fullfile(out,'figure');end
if ~isfolder(out),mkdir(out);end
% Digitized coordinates are embedded below; source PDFs are not required.
saved_path=path;path_cleanup=onCleanup(@()path(saved_path)); %#ok<NASGU>
addpath(fullfile(root,'DEA_full','baseline','code'),'-begin');
if ~isfolder(figure_out),mkdir(figure_out);end

% ====== Digitized marker centres and axis calibration ======
x_pdf=[361.7535;374.4755;387.3895;400.1115;412.9295;426.4165; ...
       438.5655;452.0525;464.2005;496.8185;528.1935];
y_pdf=[719.3545;716.1985;711.1280;704.8145;698.5980;690.3705; ...
       682.0490;673.7280;663.3965;634.3180;599.1165];
x_ticks=[367.486 399.546 431.531 463.514 495.574 527.559];
y_ticks=[721.438 683.935 646.431 608.852 571.349];
U_calibration=polyfit(x_ticks,25:5:50,1);
I_calibration=polyfit(y_ticks,0:5:20,1);
U_kV=polyval(U_calibration,x_pdf);
I_uA_cm=polyval(I_calibration,y_pdf);
U_fig5_kV=[24.107044986712747;26.09390797347781;28.098279643167402; ...
 30.102651312856988;32.10702298254658;34.19856554169019;36.09835077029338; ...
 38.18980019814485;40.08967855804019;45.196532963180445;50.11163003778029];
I_fig5_uA_cm=[.2762880689616103;.6914914710843029;1.3823042989837069; ...
 2.2141264063391004;3.0379098569535654;4.1457773094577925;5.24883271654378; ...
 6.352881392848173;7.738274359300293;11.603232886628629;16.299900370577404];
assert(numel(U_kV)==11 && all(diff(U_kV)>0) && all(I_uA_cm>0));
points=table(U_kV,I_uA_cm,x_pdf,y_pdf,U_fig5_kV,I_fig5_uA_cm);

% ====== Fixed n=2 and onset fit, plus sensitivity alternatives ======
n=2;U0_kV=21.6; % Negative-corona empirical onset in the cited paper.
x=(U_kV-U0_kV).^n;
K=sum(x.*I_uA_cm)/sum(x.^2);
fits=repmat(metrics('n2_fixed_onset_current_LS',K,n,U0_kV,U_kV,I_uA_cm),5,1);
sqrt_fit=polyfit(U_kV,sqrt(I_uA_cm),1);
fits(2)=metrics('n2_sqrt_linear_LS',sqrt_fit(1)^2,2,-sqrt_fit(2)/sqrt_fit(1),U_kV,I_uA_cm);
[u0_free,k_free]=free_onset(U_kV,I_uA_cm);
fits(3)=metrics('n2_free_onset_current_LS',k_free,2,u0_free,U_kV,I_uA_cm);
log_fit=polyfit(log10(U_kV-U0_kV),log10(I_uA_cm),1);
fits(4)=metrics('free_n_fixed_onset_log_LS',10^log_fit(2),log_fit(1),U0_kV,U_kV,I_uA_cm);
K_log=exp(mean(log(I_uA_cm)-2*log(U_kV-U0_kV)));
fits(5)=metrics('n2_fixed_onset_log_LS',K_log,2,U0_kV,U_kV,I_uA_cm);

% ====== Current-law units and two distinct geometry-transfer assumptions ======
p=read_model_params(root);
a_source_m=.0007;b_source_m=.1035;
R0_m=p.r0_m+p.d0_m;
B_source_m=a_source_m*log(b_source_m/a_source_m);
B_current_m=R0_m*log(p.rg_m/R0_m);
for k=1:numel(fits)
 exponent=fits(k).n;
 % 1 microampere/cm=1e-4 A/m; 1 kV=1000 V.
 K_SI=fits(k).K_uA_cm_kVn*1e-4/1000^exponent;
 fits(k).K_SI_A_per_m_Vn=K_SI;
 % Baseline: transfer the field coefficient inferred in the source geometry.
 fits(k).C_source_field=K_SI*B_source_m^exponent/(2*pi*a_source_m);
 % Alternative: transfer the voltage coefficient instead. Not the baseline C.
 fits(k).C_current_from_voltage_law=K_SI*B_current_m^exponent/(2*pi*R0_m);
 fits(k).I_initial_voltage_transfer_uA_m=K_SI*max(abs(p.U2_V)-p.Eth_V_m*B_current_m,0)^exponent*1e6;
 fits(k).I_initial_field_transfer_uA_m=fits(k).C_source_field* ...
  max(abs(p.U2_V)/B_current_m-p.Eth_V_m,0)^exponent*2*pi*R0_m*1e6;
end
fit_table=struct2table(fits);

% ====== Verify the selected C, n and conductivity in actual params.m ======
C_fitted=fits(1).C_source_field;
C_rounded=round(C_fitted/1e-15)*1e-15;
kappa_reference_S_m=3.23e-14;E_B_V_m=40e6;
assert(p.current_n==2,'This verification targets the selected n=2 baseline.');
assert(abs(p.spring.current_C-C_rounded)<1e-27 && abs(p.gas.current_C-C_rounded)<1e-27, ...
 'Model C does not match the selected rounded field-transfer fit.');
assert(p.spring.kappa_S_m==kappa_reference_S_m && p.gas.kappa_S_m==kappa_reference_S_m, ...
 'Model conductivity differs from the selected literature reference.');
model_check=struct('n',p.current_n,'C_fitted_A_V2',C_fitted,'C_selected_A_V2',p.spring.current_C, ...
 'C_rounding_difference_pct',100*(p.spring.current_C/C_fitted-1), ...
 'kappa_reference_S_m',kappa_reference_S_m,'kappa_inputs_match',true);
maximum_field=0;
folders={fullfile(root,'DEA_full','baseline','results'), ...
 fullfile(root,'DEA_full','sensitive','results','spring_radius')};
for j=1:numel(folders)
 saved=load_results(folders{j});
 for k=1:numel(saved.cases)
  T=saved.cases{k};maximum_field=max(maximum_field,max(abs(T.V_DEA_V)./T.thickness_m));
 end
end
model_check.max_membrane_field_MV_m=maximum_field/1e6;
model_check.max_conductivity_exponential_factor=exp(maximum_field/E_B_V_m);
model_check.ohmic_limit_note='Constant kappa is a low-field approximation; this factor is not a bound on output error.';
model_check.initial_corona_uA_m=p.spring.current_C* ...
 max(abs(p.U2_V)/B_current_m-p.Eth_V_m,0)^2*2*pi*R0_m*1e6;
result=struct('source_doi','10.1063/1.4907234','source_current_unit','microampere per centimetre', ...
 'digitization_note','Vector marker centres previously extracted from Figs. 3 and 5; not raw measurements or automatic MATLAB PDF extraction.', ...
 'B_source_m',B_source_m,'B_current_m',B_current_m, ...
 'max_U_fig3_fig5_difference_kV',max(abs(U_kV-U_fig5_kV)), ...
 'max_I_fig3_fig5_difference_uA_cm',max(abs(I_uA_cm-I_fig5_uA_cm)), ...
 'fits',fits,'model_check',model_check);
writetable(points,fullfile(out,'negative_corona_points.csv'));
writetable(fit_table,fullfile(out,'fit_comparison.csv'));
save(fullfile(out,'verification.mat'),'result','points','fit_table');
fid=fopen(fullfile(out,'verification.json'),'w');assert(fid>=0);
fprintf(fid,'%s',jsonencode(result));fclose(fid);
plot_fit(U_kV,I_uA_cm,fits,figure_out);
disp(fit_table(:,{'name','n','U0_kV','K_uA_cm_kVn','C_source_field','R2_current'}));
fprintf('Selected baseline: n=2, K=%.12g microampere/(cm kV^2)\n',K);
fprintf('Field-transfer C_fit=%.12g A/V^2; model C=%.12g A/V^2 (rounding difference %.3f%%)\n', ...
 C_fitted,p.spring.current_C,model_check.C_rounding_difference_pct);
fprintf('Both support conductivities match %.12g S/m (literature reference, not specimen validation).\n',kappa_reference_S_m);
fprintf('PASS: fit, unit/geometry conversion and selected model inputs. Output: %s\n',out);
end

function f=metrics(name,K,n,U0,U,I)
% Calculates fit errors in the original current coordinate, not log-space COD.
prediction=K*max(U-U0,0).^n;error_value=prediction-I;
f=struct('name',name,'K_uA_cm_kVn',K,'n',n,'U0_kV',U0, ...
 'rmse_uA_cm',sqrt(mean(error_value.^2)), ...
 'R2_current',1-sum(error_value.^2)/sum((I-mean(I)).^2), ...
 'mape_pct',100*mean(abs(error_value./I)), ...
 'max_abs_relative_pct',100*max(abs(error_value./I)));
end

function [U0,K]=free_onset(U,I)
% Profile K analytically and search the remaining onset coordinate.
lo=0;hi=min(U)-1e-6;g=(sqrt(5)-1)/2;
for k=1:150
 c=hi-g*(hi-lo);e=lo+g*(hi-lo);
 if loss(c,U,I)<loss(e,U,I),hi=e;else,lo=c;end
end
U0=(lo+hi)/2;x=(U-U0).^2;K=sum(x.*I)/sum(x.^2);
end

function value=loss(U0,U,I)
x=(U-U0).^2;K=sum(x.*I)/sum(x.^2);value=sum((K*x-I).^2);
end

function p=read_model_params(root)
% Resolve this package's params.m even if another model is on MATLAB's path.
saved_folder=pwd;saved_path=path;
cleanup=onCleanup(@()restore_location(saved_folder,saved_path)); %#ok<NASGU>
cd(fullfile(root,'DEA_full','baseline','code'));clear params;
p=params();
end

function restore_location(folder,saved_path)
cd(folder);path(saved_path);
end

