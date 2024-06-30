function psidata = PreprocessPsiData_dc()
%% read data from Psi-tri output files
% all the points or vertexes coordinates
data = load('D:\Psi-Tri\Psi-Tri\IDCD-run\Psitri_dci.mat');

psidata.npsi = double(data.npsi);
psidata.ntheta = double(data.ntheta);%get only z>0 half data
psidata.R = flip(data.R); 
psidata.Z = flip(data.Z); 
% the psi value of all the points/vertexes
psidata.psi = flip(data.psi'*ones(1,data.ntheta+1));

%% preprocess of the raw data
% renormalize the psi data in order to consist with ERATO framework
psi_p = max(psidata.psi(:));
psidata.psi = (psi_p-psidata.psi)/psi_p;
