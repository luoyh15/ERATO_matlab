function [TR,psidata] = PreprocessPsiData()
%% read data from Psi-tri output files
% all the points or vertexes coordinates
points = h5read('D:\Psi-Tri\Psi-Tri\IDCD-run\mesh.h5','/R');
points = points';
% the three vertexes numbering of all the triangulars
triangulars = h5read('D:\Psi-Tri\Psi-Tri\IDCD-run\mesh.h5','/LC');
triangulars = triangulars'+1;% note that matlab's numbering start from 1
% the psi value of all the points/vertexes
psidata = h5read('D:\Psi-Tri\Psi-Tri\IDCD-run\scalar_dump.h5','/Psi0018');

%% preprocess of the raw data
% construct the triangulation mesh using matlab's triangulation structure
TR = triangulation(double(triangulars),points);
% number of all the triangulars
n_triangulars = size(triangulars,1);
% renormalize the psi data in order to consist with ERATO framework
psi_p = max(psidata);
psidata = (psi_p-psidata)/psi_p;