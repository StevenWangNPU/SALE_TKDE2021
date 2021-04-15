% Given the distance to construct a sparse graph
function  S =  construct_S( dis, idx, h, r, ni)
S = zeros(ni,ni);
a = zeros(ni,ni);
for j = 1:ni
    a(j) = sum(dis(j,2:h+1).^(1/(1-r)));
    for k = 1:h
        S(j, idx(j,k+1)) = (dis(j,k+1)).^(1/(1-r))/a(j);
    end
end
