srange = (1:6).^2*10;
chirange = srange+1;
% gr = zeros(length(srange),length(chirange));
gr = zeros(length(srange));
mu0 = -1.5;
for i=1:length(srange)
%     for j=1:length(chirange)
        n_s = srange(i);
        n_chi = chirange(i);
        GetGrowthrate;
        gr(i) = -GrowthRate;
        mu0 = GrowthRate;
%     end
end
% [sr,chir]=meshgrid(srange,chirange);
% figure(2);mesh(sr,chir,-gr);
figure(3);
plot(srange,gr);
xlabel('$$N_s$$','Interpreter','latex');
ylabel('$$\Gamma^2$$','Interpreter','latex');