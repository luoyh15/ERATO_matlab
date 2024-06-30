function q = Getqprofiles()
%% read data from Psi-tri output files
data = readtable('safety_factor.dat');
data = table2array(data);
q.psi = data(:,1);
% q.psi_p = max(psi);
% q.psi = (q.psi_p-q.psi)./q.psi_p;
q.q = data(:,2);
