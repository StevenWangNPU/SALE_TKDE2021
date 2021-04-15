function  S0 =  construct_S4( idx, h,n)
% initial a graph in which the weight of nearest
% neighborhood is one.
S0 = zeros(n,n);
for j = 1:n
    for k = 1:h
        S0(j, idx(j,k+1)) = 1/h;
    end
end