function [psi_p,fT_psi,fp_psi,fq_psi] = GetGSprofiles()
%% read data from Psi-tri output files
% fileID_gs = fopen('gs.prof','r');
% data = fscanf(fileID_gs,'%f',[5,Inf]);
% fclose(fileID_gs);
% psi = data(1,:)';
% gs.psi_p = max(psi);
% gs.psi = (gs.psi_p-psi)./gs.psi_p;
% gs.dTDpsi = -data(2,:)';
% gs.T = data(3,:)';
% gs.dPDpsi = -data(4,:)';
% gs.P = data(5,:)';
data = load('D:\Psi-Tri\Psi-Tri\IDCD-run\Psitri_dci.mat');
% psi_p = max(data.psi);
psi_p = 1;
psi = (psi_p-data.psi)/psi_p;
fT_psi =  fit(psi',data.T','spline');
% T0 = fT_psi(0);
% fT_psi = fit(psi',data.T'/T0,'spline');
fp_psi = fit(psi',data.P','spline');
fq_psi = fit(psi',data.q','spline');

