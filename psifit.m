p_all = [];
psi_all = [];
for i=1:length(pF)
    for j = 2:length(pF{i}.functions)
        f=pF{i}.functions(j-1,:);
        p=pF{i}.points(j,:);       
        psi_all = [psi_all;f*[1,p(1),p(2),p(1).^2,p(1).*p(2),p(2).^2]'.^2];
        p_all = [p_all;p];
    end
end
