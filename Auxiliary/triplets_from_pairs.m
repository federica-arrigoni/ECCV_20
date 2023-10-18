
function [triplets,n_triplets]=triplets_from_pairs(A,c)
% A = viewing graph
% ncams = number of cameras
% c = number of triplets sampled at random from each pair
% triplets = cell containing sampled triplets
% n_triplets = number of sampled triplets

% Find all the possible image pairs
[I,J]=find(triu(A,1));

% for each pair, c triplets are sampled at random
pairs=repmat([I J],c,1);
n_pairs=size(pairs,1);

triplets=[]; % initialize
for k=1:n_pairs
    
    i=pairs(k,1); j=pairs(k,2);
    
    % find the edges incident to node i
    nodes_i=find(A(i,:));
    % find the edges incident to node j
    nodes_j=find(A(j,:));
    
    % edges incident to edge (i,j)
    candidates=[nodes_i nodes_j];
    candidates(candidates==i)=[];
    candidates(candidates==j)=[];
    candidates=unique(candidates); % remove duplicates
    
    ncandidates=length(candidates);
    
    if ncandidates>0
        l=randi(ncandidates);
        triplets=[triplets; sort([i j candidates(l)])];
    end
        
%     u=randi(ncams);
%     while u==triplets(k,1) || u==triplets(k,2)
%         u=randi(ncams);
%     end
%     
%     triplets(k,3)=u;
%     triplets(k,:)=sort(triplets(k,:));
end

% remove duplicates
triplets=unique(triplets,'rows');
n_triplets=size(triplets,1);

triplets=mat2cell(triplets,ones(1,n_triplets),3);

end

