function [ps,pchi,ms,mchi] = NonEqualMesh(psi_p,fq_psi,n_s,n_chi,n)
% the equal mesh 
ms_eql = (0:1/n_s:1)';
ms = zeros(size(ms_eql));% store the new s mesh
% function q of s
fq_s = @(s) fq_psi(psi_p*s.^2);
fqDs = @(s) 2*s.*differentiate(fq_psi,psi_p*s.^2);
q0 = fq_s(0);qs = fq_s(1);
m1 = ceil(n*qs); m2 = floor(n*q0);
n_sing = 0;
x_sing = [];
delta_coeff = 1;% more dense around sigular points if this is smaller
for mi = m1:m2
    xi = fzero(@(ss) n*fq_s(ss)-mi,[0,1]);
    n_sing = n_sing+length(xi);
    x_sing = [x_sing;xi];
end
if ~isempty(x_sing)
    delta_sing = delta_coeff./(n*fqDs(x_sing));
    alpha_sing = 2*n_sing/(3*(n_sing+1))...%constant to prevent too many points around sigularity
        /sum(atan((1-x_sing)./delta_sing)+atan(x_sing./delta_sing),1);
    W_sing = @(s) (3+n_sing)/(3*(n_sing+1))*s+...
        alpha_sing*sum(atan((s-x_sing)./delta_sing)+atan(x_sing./delta_sing),1);
    for i = 2:length(ms_eql)-1
        ms(i) = fzero(@(s) W_sing(s)-ms_eql(i),[0,1]);
    end
    ms(end) = 1;
else
    msg = 'no sigularity found.';
    disp(msg);
    ms = ms_eql;
end

ps = (ms(1:end-1)+ms(2:end))/2;
% chi direction is equal mesh
mchi = (-pi/(2*(n_chi-1)):pi/(n_chi-1):pi+pi/(2*(n_chi-1)))';
pchi = (0:pi/(n_chi-1):pi)';
