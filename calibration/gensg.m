function sg = gensg(n,d,varargin)
% GENSG Create a Sparse Grid Design
%
% sg=gensg(eta,d) returns a structure that constains information about a d
% demensional sparse grid design with level of constrcution eta.  The
% component designs are prechosen in an ad-hoc manor.
% sg = gensg(eta,d,X) does the same operation using component designs in
% the structure of different sized vectors X of size d x k (k>eta-d+1) where X{i+1,j} must contain
% X{i,j}.
%  This code can take awhile for large eta and/or d, please be patient.
%
%
% This code is to be used with the Sparse Grid Designs package written with
% the paper "Fast prediction of deterministic functions using sparse grid
% experimental designs" by Matthew Plumlee.
%
% Copyright, 2014, Matthew Plumlee.  All rights reserved
if nargin == 3;
    X = varargin{1}';
    X_temp = X;
    for k = 1:d
        for j = 1:length(X)
            if j >1
                X{j,k} = setdiff(X_temp{j,k},X_temp{j-1,k});
            end
            X_tot{j,k} = sort([X{1:j,k}]);
            X_dim(j,k) = length(X{j,k});
            X_dim_tot(j,k) = length(X_tot{j,k});
        end
    end
else
    for k = 1:d
        X{1,k} =   .5   ;
        %X{2,k} = [ .125,.875 ];
        X{2,k}=[ .25 .75 ];
        %X{3,k} = [.25,.75];
        X{3,k}=[.125 .375 .625 .875];
        %X{4,k} = [0,1];
        X{4,k}=[.0625 .1875 .3125 .4375 .5625 .6875 .8125 .9375];
        %X{5,k} = [.375,.625];
        X{5,k}= 1/32:2/32:(1-1/32);
        X{6,k} = [.25-1/16,.75+1/16];
        X{7,k} = [1-1/16,0+1/16];
        X{8,k} = [.5-1/16,.5+1/16];
        X{9,k} = [0.3125,   0.6875];
        X{10,k} = [.25+1/32,.75-1/32];
        X{11,k} = [0+1/32,1-1/32];
        X{12,k} = [.875+1/32,.125-1/32];
        X{13,k} = [.875-1/32,.125+1/32];
        X{14,k} = [.375+1/32,.625-1/32];
        X{15,k} = [.5-1/32,.5+1/32];
        X{16,k} = [.25-1/32,.75+1/32];
        for j = 1:length(X)
            X_dim(j,k) = length(X{j});
            X_tot{j,k} = sort([X{1:j}]);
            X_dim_tot(j,k) = length(X_tot{j});
        end
    end
end
totlevelsets = 0;
numlevelset = zeros(n-d,2);
for i = d:n
    totlevelsets = totlevelsets+nchoosek(i-1,d-1);
    numlevelset(i-d+1,:) = [totlevelsets-nchoosek(i-1,d-1)+1 totlevelsets];
end
totlevelset =ones(totlevelsets,d);
k = 0;
for i = d:n
    totlevelset(k+1:k+nchoosek(i-1,d-1),:) = spgetseq((i-d),d)+1;
    k = k+nchoosek(i-1,d-1);
end
coeff = zeros(size(totlevelset,1),1);
for j = 1:size(totlevelset,1)
    abs_j = sum(totlevelset(j,:));
    if max(d,n-d+1)<= abs_j
        coeff(j) = (-1)^(n-abs_j)*nchoosek(d-1,n-abs_j);
    end
end
num_points = 0;
for j  = 1:size(totlevelset,1)
    q_inds= sub2ind(size(X_dim), totlevelset(j,:), 1:d);
    num_points = num_points + prod(X_dim(q_inds));
end
X_points= zeros(num_points,d);
C_mat = zeros(n,d,n-d);
for k= 1:n
    for i = 1:d
        for j = 1:(n-d)
            C_mat(k,i,j) = counting(k,i,j);
        end
    end
end
k=0;
clear idx
X_points_temp = zeros(max(X_dim(n-d+1,:)),d);
label = cell(size(totlevelset,1),1);
ulabel = cell(size(totlevelset,1),1);
ind_store2 = cell(size(totlevelset,1),1);
for j = 1:size(totlevelset,1)
    for l = 1:(n-d+1)
        if numlevelset(l,1) <= j && numlevelset(l,2) >= j
            break;
        end
    end
    q_inds= sub2ind(size(X_dim), totlevelset(j,:), 1:d);
    
    q = X_dim(q_inds);
    X_points_temp(1:prod(q),d) = repmat(X{totlevelset(j,d),d}',prod(q(1:d-1)),1);
    b = prod(q(2:d))*ones(1,q(1));
    idx([cumsum([1 b(b>0)])]) = 1;
    X_points_temp(1:prod(q),1) = X{totlevelset(j,1),1}(cumsum(idx(1:find(idx,1,'last')-1)))';
    clear idx
    for i = 2:d-1
        b = prod(q(i+1:d))*ones(1,q(i));
        idx([cumsum([1 b(b>0)])]) = 1;
        X_points_temp(1:prod(q),i) = repmat(X{totlevelset(j,i),i}(cumsum(idx(1:find(idx,1,'last')-1))),1,prod(q(1:i-1)));
        clear idx
    end
    X_points(k+1:k+prod(q),:) = X_points_temp(1:prod(q),:);
    
    f = X_dim_tot(q_inds);
    
    label{j} = zeros(prod(f),1);
    k_0 = 0;
    ulabel{j} = (k+1):(k+prod(q));
    
    if l>1
        [~,vals] = find(totlevelset(j,:)>=2);
        indexed_vals =zeros(1,length(vals));
        if l >2
            clear indexed_vals
            for inter = 1:length(vals)
                samp_level_set = totlevelset(j,:);
                samp_level_set(vals(inter)) = totlevelset(j,vals(inter))-1;
                indexed_vals(inter) = counter_f(samp_level_set,C_mat)+numlevelset(l-1,1)-1;
            end
            ind_store2{j} = sort([unique([ind_store2{indexed_vals}]),indexed_vals]);
        else
            indexed_vals = 1;
            ind_store2{j} = indexed_vals;
        end
        
        for j_0 = ind_store2{j}
            label{j}(k_0+1:k_0+length(ulabel{j_0})) = ulabel{j_0};
            k_0 = k_0 + length(ulabel{j_0});
        end
    end
    label{j}(k_0+1:end) = ulabel{j};
    k = k+prod(q);
    [~,ind_val] = sortrows(X_points(label{j},:),1:d);
    label{j} = label{j}(ind_val);
end
sg.design = X_points;
sg.labels = label;
sg.grids = X_tot;
sg.eta = n;
sg.d = d;
sg.coeff = coeff;
sg.levelsets = totlevelset;
sg.griddim = X_dim_tot;
end
function index_vals_temp = counting(n,d,c)
index_vals_temp=0;
warning off
if d-2>=0
    if n-c>=d-2
        for lcv = 2:1:c
            index_vals_temp =  index_vals_temp+nchoosek(double(n-lcv),double(d-2));
        end
    end
end
warning on
end
function k = counter_f(vec_val,C_mat)
k=1;
n= sum(vec_val);
d= length(vec_val);
for i = d:-1:2
    if vec_val(i)>1
        k = k+ C_mat(n-sum(vec_val(i+1:end)),i,vec_val(i));
    end
end
end
