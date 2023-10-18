
function [tracks_triplets,X_triplets,labels_triplets,triplets,n_triplets] = triplet_segmentation(pairwiseEst,SIFT,method_tracks,triplets,n_triplets,d,sigma_gt)

if nargin<7
    sigma_gt=0.1;
end

%% Compute tracks

tracks_triplets=cell(n_triplets,1);
X_triplets=cell(n_triplets,1);
labels_triplets=cell(n_triplets,1);

parfor tt=1:n_triplets % parallel computation to speed up
    
    [tt n_triplets]
    
    %% Tracks
    
    tri=triplets{tt};
    [XX,TT] = tracks_3view(pairwiseEst,tri(1),tri(2),tri(3),SIFT,method_tracks);
    size(XX)
    
    %%
    if size(XX,2)>15
        %% Three-frame segmentation
        LL=modelfitting_3view(XX(1:2,:),XX(4:5,:),XX(7:8,:),sigma_gt,d);
        
        %% Save data
        
        X_triplets{tt}=XX;
        tracks_triplets{tt}=TT;
        labels_triplets{tt}=LL;
        
    else
        X_triplets{tt}=[];
        tracks_triplets{tt}=[];
        labels_triplets{tt}=[];
    end
    
end

%% Remove empty (failure) three-frame segmentations

ind_empty=[];
for tt=1:n_triplets
    if isempty(labels_triplets{tt})
        ind_empty=[ind_empty tt];
    end
end

labels_triplets(ind_empty)=[];
tracks_triplets(ind_empty)=[];
X_triplets(ind_empty)=[];
triplets(ind_empty)=[];
n_triplets=length(triplets);


end