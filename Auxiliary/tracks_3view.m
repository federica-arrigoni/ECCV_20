
function [X,tracks] = tracks_3view(pairwiseEst,i,j,k,SIFT,method)

% compute matches
if strcmp(method,'consistent')
    tracks = tracks_triplet(pairwiseEst,i,j,k,method);
else
    tracks1 = tracks_triplet(pairwiseEst,i,j,k,method);
    tracks2 = tracks_triplet(pairwiseEst,j,k,i,method); tracks2=[tracks2(:,3) tracks2(:,1) tracks2(:,2)];
    tracks3 = tracks_triplet(pairwiseEst,k,i,j,method); tracks3=[tracks3(:,2) tracks3(:,3) tracks3(:,1)];
    tracks=unique([tracks1;tracks2;tracks3],'rows');
end

% compute points
if ~isempty(tracks)
    ind1=tracks(:,1); ind2=tracks(:,2); ind3=tracks(:,3);
    
    npoints=size(tracks,1);
    X=ones(9,npoints);
    
    SIFT_i=SIFT{i};
    Xi=[SIFT_i.locs(ind1,1)';SIFT_i.locs(ind1,2)'];
    X(1:2,:)=Xi;
    
    SIFT_j=SIFT{j};
    Xj=[SIFT_j.locs(ind2,1)';SIFT_j.locs(ind2,2)'];
    X(4:5,:)=Xj;
    
    SIFT_k=SIFT{k};
    Xk=[SIFT_k.locs(ind3,1)';SIFT_k.locs(ind3,2)'];
    X(7:8,:)=Xk;
else
    X=[];
end

end


function tracks = tracks_triplet(pairwiseEst,i,j,k,method)

tracks=[];

if i<j
    ind_i_ij=pairwiseEst{i,j}.ind1;
    ind_j_ij=pairwiseEst{i,j}.ind2;
else
    ind_i_ij=pairwiseEst{j,i}.ind2;
    ind_j_ij=pairwiseEst{j,i}.ind1;
end

if j<k
    ind_j_jk=pairwiseEst{j,k}.ind1;
    ind_k_jk=pairwiseEst{j,k}.ind2;
else
    ind_j_jk=pairwiseEst{k,j}.ind2;
    ind_k_jk=pairwiseEst{k,j}.ind1;
end

if i<k
    ind_i_ik=pairwiseEst{i,k}.ind1;
    ind_k_ik=pairwiseEst{i,k}.ind2;
else
    ind_i_ik=pairwiseEst{k,i}.ind2;
    ind_k_ik=pairwiseEst{k,i}.ind1;
end


for p=1:length(ind_i_ij) % matches in the pair (i,j)
    
    ind_i=ind_i_ij(p); % point in image i
    ind_j=ind_j_ij(p); % point in image j
    
    [a,b]=ismember(ind_j,ind_j_jk);
    
    if a % if the point is present also in image k...
        
        ind_k=ind_k_jk(b); % point in image k
        
        % check for consistency...
        [c,d]=ismember(ind_i,ind_i_ik);
        [e,f]=ismember(ind_k,ind_k_ik);
        
        if strcmp(method,'consistent') % keep only consistent tracks
            
            if c && e % consistency
                tracks=[tracks;ind_i ind_j ind_k];
            end
        
        elseif strcmp(method,'remove_inconsistent') % discard inconsistent tracks
            
            if c && e % consistency
                tracks=[tracks;ind_i ind_j ind_k];
                disp('the track is CONSISTENT')
            elseif c || e % inconsistency
                disp('the track is not consistent')
            else % the point is not matched between images i and k
                tracks=[tracks;ind_i ind_j ind_k];
                disp('the track is not checked')
            end
            
        elseif strcmp(method,'all') % keep all tracks
            tracks=[tracks;ind_i ind_j ind_k];
        else
            warning('unknown method')
        end
        
    end
end

end
