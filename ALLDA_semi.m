function [W,p,S,Obj]  = ALLDA_semi( LX,Y,X,h1,h2,m,alpha,maxiter)
% LX: labeled data, every column is a sample
% Y:  label vector
% X:  all train data  
% h1: h1 nearest neighborhood in labeled data
% h2: h2 nearest neighborhood in all data
% m:  final dimensions
% alpha: a parameter
% maxiter: the maximal numbre of iteration
% written by Steven Wang
% email: zhengwangml@gmail.com

%% load scale
[~,n] = size(X);  % n denotes the number of all samples
c = unique(Y);    
l = length(Y);    % l denotes the number of labeled samples
p0 = [];

%% initialization
for i=1:length(c)
    Xc{i} = LX(:,Y == i);  %  store the each class samples
    nc(i) = size(Xc{i},2);  %  record the number of each class
end
H = eye(l) - 1/l*ones(l);
% H = eye(80) - 1/80*ones(80);
St = LX'*H*LX;  % St of labeled data
invSt = inv(St);
% initialized P similarity matrix of labeled data
for k = 1 : length(c)
    Xi = Xc{k};
    ni = nc(k);   
    distXi = L2_distance_1(Xi,Xi);
    [~, idx] = sort(distXi,2);
    S0{k} = construct_S2(idx, h1, ni);   
    p0 = blkdiag(p0,S0{k});         
end
obj1 = zeros(1,length(c));
Obj = zeros(1,maxiter);

% initialized S similarity matrix of all data
distXX = L2_distance_1(X,X);
[~, idx] = sort(distXX,2);
S0 = construct_S4(idx,h2,n);

%% Iteration
for iter = 1:maxiter
p = [];
P = p0;
S = S0.^2;

% Calculate lapalcian matrix L_p;
P = (P+P')/2;
D_p = diag(sum(P));
L_p = D_p - P;

% Calculate lapalcian matrix L_s;
S = (S+S')/2;
D_s = diag(sum(S));
L_s = D_s - S;

% Calculate projection matrix W
G = invSt*(LX*L_p*LX'+alpha*X*L_s*X'); 
[W,~,~] = eig1(G, m, 0, 0);
W = W*diag(1./sqrt(diag(W'*St*W)));

% Updata matrix P
 for i = 1:length(c)
 Xc{i} = LX(:,Y==i); 
 nc(i) = size(Xc{i},2);
 Xi = Xc{i};
 ni = nc(i);
 distLXx = L2_distance_1(W'*Xi,W'*Xi);
 [~, idx] = sort(distLXx,2);
 PP{i} = construct_S2( idx, h1, ni);
 p = blkdiag(p,PP{i});
 obj1(i) = sum(sum(PP{i}.*distLXx));
% obj2(i) = sum((sum (dis(:,2:h+1).^(1/(1-2)),2) ).^(1-2));
 end
Obj1 = sum(obj1); 

% Updata matrix S
distXx = L2_distance_1(W'*X,W'*X);
[disX, idxX] = sort(distXx,2);
S =  construct_S(disX + eps, idxX, h2, 2, n);
Obj2 = alpha*sum(sum(distXx.*(S.^2)));

% Change the variable
 p0 = p;
 S0 = S;
 
%  Objective function value
Obj(1,iter) = Obj1 + Obj2;

end


