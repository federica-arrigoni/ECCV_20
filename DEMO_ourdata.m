
clear,clc,close all
addpath('./Auxiliary/')
addpath(genpath('./RPA_beta/'))

%% Load data

folder_path = './DATASETS_ECCV20/';
dataset='STUFFED_ANIMALS6';

load([folder_path dataset '/RESULTS.mat'])
load([folder_path dataset '/labels_gt.mat'])
method_permutation='hungarian';
[d ncams m]

%% TRISEG: Triplet-based segmentation

disp(['TRISEG:'])
tic
group2=segment_mode_triplets(labels_triplets,triplets,tracks_triplets,dim,ncams,d,method_permutation);
toc

% Compute error
[missrate2,known2,group2,indknown2]=compute_missrate(group2,labels_gt);
disp(['Missclassification error (classified points): ' num2str(missrate2*100) '%'])
disp(['Percentage of classified points: ' num2str(known2*100) '%'])
disp(' ')

%% TRIPAIRSEG: Consider both pairs and triplets

disp(['TRIPAIRSEG:'])
% Define adjacency matrix
for i=1:ncams
    for j=i+1:ncams
        if isempty(labels_pairwise{i,j})
            A(i,j)=0; A(j,i)=0;
        end
    end
end

[I,J]=find(triu(A,1));
npairs=length(I);
pairs=mat2cell([I J],ones(1,npairs),2);

labels_pairs=cell(npairs,1);
tracks_pairs=cell(npairs,1);

for kk=1:npairs
    i=I(kk); j=J(kk);
    labels_pairs{kk}=labels_pairwise{i,j};
    tracks_pairs{kk}=[pairwiseEst{i,j}.ind1' pairwiseEst{i,j}.ind2'];
end
subsets=[pairs;triplets];
labels_subsets=[labels_pairs;labels_triplets];
tracks_subsets=[tracks_pairs;tracks_triplets];

tic
group3=segment_mode_general(labels_subsets,subsets,tracks_subsets,dim,ncams,d,method_permutation);
toc

% Compute error
[missrate3,known3,group3,indknown3]=compute_missrate(group3,labels_gt);
disp(['Missclassification error (classified points): ' num2str(missrate3*100) '%'])
disp(['Percentage of classified points: ' num2str(known3*100) '%'])
disp(' ')


