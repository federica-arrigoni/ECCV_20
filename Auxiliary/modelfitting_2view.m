

function group=modelfitting_2view(y1,y2,ngroups,method)

if nargin<4
    method='RPA';
    %method='RansaCov';
end

%% Parameters

sigma_gt=0.005; % ICCV, BMVC


%% Set model

model = 'fundamental';
%model = 'homography';
[ distFun, hpFun, fit_model, cardmss] = set_model( model );

%% normalization

npoints=size(y1,1);
y1=[y1';ones(1,npoints)];
y2=[y2';ones(1,npoints)];

[dat_img_1 T1] = normalise2dpts(y1);
[dat_img_2 T2] = normalise2dpts(y2);
X = [ dat_img_1 ; dat_img_2 ];

%% sampling hypotheses

% Ransacov
% S = mssWeighted( X, 3*npoints, 1*npoints,'loc', model, 0.6, nan );

% guided sampling
w = 0.5;
blk = 3*npoints;
S  = mssWeighted( X, 6*npoints, blk, 'cauchy', model, w, sigma_gt);
S=S(blk+1:end,:);

H = hpFun(X,S); %hypotheses

R = res( X, H, distFun ); disp('Residuals computed')

%%

if strcmp(method,'RPA')
    
    group=model_fittingRPA(X,S,R,model,sigma_gt,ngroups);
    
elseif strcmp(method,'Tlinkage')
    
    disp('....')
    
    % Preference Matrix
    P = prefMat(R, sigma_gt, 6);
    
    % Clustering
    group = tlnk(P);
    
    % Outlier rejection step
    group  = outlier_rejection_card( group, cardmss );
    
    disp([num2str( max(group) ),' motions detected'])
    
elseif strcmp(method,'RansaCov')
    
    epsilon = .4e-1; % inlier threshold
    
    P = zeros(size(R));
    P(R<epsilon) = 1;
    
    % preprocessing step (optional)
    [Q, ~] = refit_consensus( X, P, H, epsilon, model );
    F = pruning_consensus(Q);
    
    %select the k structures covering more points
    ind_group = ransacov(F,ngroups);
    
    ind_group
    
    group=zeros(npoints,1);
    
    ind_inters_out=find(sum(ind_group')~=1); % intersecting structures or outliers
    ind_group(ind_inters_out,:)=0;
    
    ind_group
    
    [I,J]=find(ind_group);
    group(I)=J;
end


end


function group=model_fittingRPA(X,S,R,model,sigma_gt,ngroups)


%% Preference Trick
P = prefMat(R, sigma_gt, 6); % preference matrix

K = exp(- (squareform(pdist(P,@tanimoto))).^2);  % similarity matrix

%% Robust PCA
try
    lambda = 1/sqrt(size(K,1));
    [K_rpca, E_hat, ~] = inexact_alm_rpca(K, lambda);
    
    
    %% symmetric matrix factorization
    
    %[Uinit, mekmeans]  = guess_init( K_rpca, ngroups , G);
    [Uinit]  = guess_init( K_rpca, ngroups);
    
    params.Hinit=Uinit; params.maxiter = 100000;
    [U, iter, obj] = symnmf_anls(K_rpca, ngroups,  params);
    indU = indMax(U);
    
    % segmentations obtained from snmf
    % NB: F is a segmentation
    F = seg_from_binaryU(U);
    softIndU= U;    softIndU(indU==0)=0; % mlsac like
    
    %try
    %% Model extraction
    
    niter = 1000;
    [Phi, Z] =  rinforzino(X, S, P, F, softIndU , model, sigma_gt, niter);
    
    [~,I] = max(Phi'*softIndU,[],1);
    mss = Z(I,:);
    
    
    %% refinement using robust statistic
    
    cost = 1.5271;
    group = segmentation( mss, X, model, U , sigma_gt, cost ,'nearest');
    
catch
    
    disp('Error in RPA: segmentation not performed')
    group=[];
    
end


end