

function group=modelfitting_3view(y1,y2,y3,sigma_gt,ngroups)


%% Parameters
 
%sigma_gt=0.1; 

%% normalization - togliere???

npoints=size(y1,2);
y1=[y1;ones(1,npoints)];
y2=[y2;ones(1,npoints)];
y3=[y3;ones(1,npoints)];

%X=[y1;y2;y3];

[dat_img_1 T1] = normalise2dpts(y1);
[dat_img_2 T2] = normalise2dpts(y2);
[dat_img_3 T3] = normalise2dpts(y3);
X = [ dat_img_1 ; dat_img_2; dat_img_3 ];

%% Set model

model = 'trifocal';
[ distFun, hpFun, fit_model, cardmss] = set_model( model );

%% sampling hypotheses

% guided sampling
w = 0.5;
blk = 3*npoints;
S  = mssWeighted( X, 6*npoints, blk, 'cauchy', model, w, sigma_gt);
S=S(blk+1:end,:);

H = hpFun(X,S); %hypotheses

R = res( X, H, distFun ); disp('Residuals computed')

%% segmentation

group=model_fittingRPA(X,S,R,model,sigma_gt,ngroups);


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

    
    %% Model extraction
    
    niter = 1000; % 10 o 1000????
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