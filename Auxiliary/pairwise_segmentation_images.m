
function [Z,labels_pairwise]=pairwise_segmentation_images(Zmatch,SIFT,dim,d,n,method,A)

if nargin<7
    A=ones(n);
end

cumDim = [0;cumsum(dim(1:end-1))];
m=sum(dim);

Z=sparse(m,m);
labels_pairwise=cell(n);

for i=1:n
    
    SIFT_i=SIFT{i};
    
    for j=i+1:n
        
        if A(i,j)==1
            
            %%
            SIFT_j=SIFT{j};
            
            [i j]
            
            Zij=Zmatch(1+cumDim(i):cumDim(i)+dim(i),1+cumDim(j):cumDim(j)+dim(j));
            
            % ind1 and ind2 encode the matches between image_i and image_j
            [ind1,ind2]=find(Zij);
            
            % xi,xj are the matching points for the pair (i,j)
            Xi=[SIFT_i.locs(ind1,1)';SIFT_i.locs(ind1,2)'];
            Xj=[SIFT_j.locs(ind2,1)';SIFT_j.locs(ind2,2)'];
            
            if length(ind1)>8 %&& ~isequal([i j],[1 5]) % if there are (at least) 8 correspondences
                length(ind1)
                if strcmp(method,'RPA') || strcmp(method,'RansaCov') || strcmp(method,'Tlinkage')
                    groups=modelfitting_2view(Xi',Xj',d,method);
                end
                
                groups'
                
                if ~isempty(groups)
                    % assign a label to ALL points in image h (zero means no label)
                    groups_i=zeros(dim(i),1);
                    groups_i(ind1)=groups;
                    
                    % assign a label to ALL points in image k (zero means no label)
                    groups_j=zeros(dim(j),1);
                    groups_j(ind2)=groups;
                    
                    % construct a binary matrix that encodes the segmentation
                    Zhk=segment2matrix(groups_i,groups_j,d);
                    
                    Z(1+cumDim(i):cumDim(i)+dim(i),1+cumDim(j):cumDim(j)+dim(j)) = Zhk;
                    Z(1+cumDim(j):cumDim(j)+dim(j),1+cumDim(i):cumDim(i)+dim(i)) = Zhk';
                    
                    % save pairwise labels
                    labels_pairwise{i,j}=groups;
                    labels_pairwise{j,i}=groups;
                end
            end
            
            
        end
    end
end


end




