
clear,clc,close all
addpath('./Auxiliary/')
addpath(genpath('./RPA_beta/'))

% Add VLFEAT to the matlab path (SIFT)
fold = pwd;
cd './vlfeat-0.9.21/toolbox/'
vl_setup
cd(fold);

% NOTE: do not add ALL VLFEAT folder to the path as it contains functions
% with the same name as in RPA that might create troubles

%% Load raw images

folder_path = './ECCV20_DATA/';
dataset='STUFFED_ANIMALS4';

img_path =[folder_path dataset '/']; format_img='jpg';

%% Parameters

%scale=1; % (no scaling)
scale=0.7; % Rescale images to speed up SIFT

imnames = dir(strcat(img_path,'*.',format_img));
ncams=length(imnames);

% delete features that have <= min_match matches
min_match=1;

d=4; % number of motions
% THIS SHOULD BE KNOWN IN ADVANCE!!!

%% Compute SIFT locations and descriptors for each image

SIFT = cell(1,ncams);
dim=zeros(ncams,1);

for i=1:ncams
    
    fprintf('\nComputing frames and descriptors: image %d \n',i);
    
    tic;
    im = imread(strcat(img_path,imnames(i).name)); % load the current image
    im = imresize(im,scale); % rescale the image to speed-up SIFT
    
    if size(im,3)==1
        im=single(im);
    else
        im=single(rgb2gray(im));
    end
    [frames1,descr1] = vl_sift(im) ; % computes SIFT locations and descriptors
    
    SIFT{i}.desc = descr1;
    SIFT{i}.locs = frames1(1:2,:)';
    SIFT{i}.locs = SIFT{i}.locs/scale;
    SIFT{i}.scale = frames1(3,:)';
    
    fprintf('%d descriptors extracted\n',size(SIFT{i}.locs,1));
    toc
    
    dim(i)=size(SIFT{i}.locs,1);
    
end

cumDim = [0;cumsum(dim(1:end-1))];


%% Match all the pairs

[~,Z_pairwise] = matching_noransac(ncams,SIFT,dim);

for i=1:ncams
    n_match=sum( Z_pairwise(1+cumDim(i):cumDim(i)+dim(i),:) ,2);
    ind_match=find(n_match<=min_match);
    
    Z_pairwise(cumDim(i)+ind_match,:)=[];
    Z_pairwise(:,cumDim(i)+ind_match)=[];
    dim(i)=dim(i)-length(ind_match);
    
    cumDim = [0;cumsum(dim(1:end-1))];
    
    SIFT{i}.desc(:,ind_match)=[];
    SIFT{i}.locs(ind_match,:)=[];
    SIFT{i}.scale(ind_match)=[];
end

m=size(Z_pairwise,1);
pairwiseEst=ZtoMatches(Z_pairwise,dim,ncams);


%% Plot matches for a selected pair

i=1; j=2;

image_i=imread(strcat(img_path,imnames(i).name)); % left image
image_j=imread(strcat(img_path,imnames(j).name)); % right image
SIFT_i=SIFT{i};
SIFT_j=SIFT{j};

Zij=Z_pairwise(1+cumDim(i):cumDim(i)+dim(i),1+cumDim(j):cumDim(j)+dim(j));
plot_sift(image_i,image_j,Zij,SIFT_i,SIFT_j);

set(gca,'FontSize',22,'LineWidth',3)
title(['Pair (' num2str(i) ',' num2str(j) ')'],'FontWeight','Normal')

%% Perform pairwise segmentation

[Z,labels_pairwise]=pairwise_segmentation_images(Z_pairwise,SIFT,dim,d,ncams,'RPA',ones(ncams)); 


%% Compute triplets

if ncams<=8 % all triplets
    triplets=nchoosek(1:ncams,3);
    n_triplets=size(triplets,1);
    triplets=mat2cell(triplets,ones(1,n_triplets),3);
else
    [triplets,n_triplets]=triplets_from_pairs(ones(ncams),2);
end

%% Perform triplet-based segmentation

% NB: all the points are used!!! tracks of length 2 are not removed
method_tracks='all';

sigma_gt=0.1; % 0.1 default
tic
[tracks_triplets,X_triplets,labels_triplets,triplets,n_triplets] = triplet_segmentation(pairwiseEst,SIFT,method_tracks,triplets,n_triplets,d,sigma_gt);
toc


%% Plot segmentation for a selected triplet

ind=3; % up to n_triplets

triplet=triplets{ind};
group_lcr=labels_triplets{ind};
tracks=tracks_triplets{ind};

l=triplet(1); c=triplet(2); r=triplet(3);

%      image_l=(imread(strcat(img_path,imnames(l).name))); % left
%      image_c=(imread(strcat(img_path,imnames(c).name))); % centre
%      image_r=(imread(strcat(img_path,imnames(r).name))); % right

image_l=rgb2gray(imread(strcat(img_path,imnames(l).name))); % left
image_c=rgb2gray(imread(strcat(img_path,imnames(c).name))); % centre
image_r=rgb2gray(imread(strcat(img_path,imnames(r).name))); % right
[rows,columns]=size(image_r);

SIFT_l=SIFT{l}; SIFT_r=SIFT{r}; SIFT_c=SIFT{c};
Xl=[SIFT_l.locs(:,1)';SIFT_l.locs(:,2)']'; % points in image l
Xr=[SIFT_r.locs(:,1)';SIFT_r.locs(:,2)']'; % points in image r
Xc=[SIFT_c.locs(:,1)';SIFT_c.locs(:,2)']'; % points in image r

ind1=tracks(:,1);
ind2=tracks(:,2);
ind3=tracks(:,3);

Xl=Xl(ind1,:); Xc=Xc(ind2,:); Xr=Xr(ind3,:);

ind_l=1+cumDim(l):cumDim(l)+dim(l); % indices of points in image l
ind_c=1+cumDim(c):cumDim(c)+dim(c); % indices of points in image c
ind_r=1+cumDim(r):cumDim(r)+dim(r); % indices of points in image r

colors = lines(7);

% Input
figure,
imshow(cat(2, image_l, image_c, image_r)) ;
hold on
set(gca,'FontSize',22,'LineWidth',3)
title(['Input - triplet (' num2str(l) ',' num2str(c) ',' num2str(r) ')'],'FontWeight','Normal')

for p=1:length(ind1)
    if group_lcr( p )~=0 % if the point has been classified
        plot(Xl(p,1),Xl(p,2),'.','Color',colors(group_lcr( p )+1,:),'MarkerSize',15)
        plot(Xc(p,1)+columns,Xc(p,2),'.','Color',colors(group_lcr( p )+1,:),'MarkerSize',15)
        plot(Xr(p,1)+2*columns,Xr(p,2),'.','Color',colors(group_lcr( p )+1,:),'MarkerSize',15)
    end
end


%% Multi-frame segmentation
%% TRIPLETS

method_permutation='hungarian';
tic
group_tri=segment_mode_triplets(labels_triplets,triplets,tracks_triplets,dim,ncams,d,method_permutation);
toc

% Classified points
ind_known=find(group_tri~=0);
known_rate=length(ind_known)/length(group_tri);
disp(['Percentage of classified points: ' num2str(known_rate*100) '%'])


%% TRIPLETS & PAIRS

A=ones(ncams);
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
group_pairs_tri=segment_mode_general(labels_subsets,subsets,tracks_subsets,dim,ncams,d,method_permutation);
toc

% Classified points
ind_known=find(group_pairs_tri~=0);
known_rate=length(ind_known)/length(group_pairs_tri);
disp(['Percentage of classified points: ' num2str(known_rate*100) '%'])


%% Plot segmentation for each image for a selected method

group_out=group_tri;
%group_out=group_pairs_tri;

colors = lines(6);

for i=1:ncams
    
    %image_i=imread(strcat(img_path,imnames(i).name)); % image i
    image_i=rgb2gray(imread(strcat(img_path,imnames(i).name))); % image i
    
    SIFT_i=SIFT{i};
    Xi=[SIFT_i.locs(:,1)';SIFT_i.locs(:,2)']'; % coordinates of points in image i
    ind_i=1+cumDim(i):cumDim(i)+dim(i); % indices of points in image i
    
    figure,
    imshow(image_i)
    hold on
    set(gca,'FontSize',22,'LineWidth',3)
    
    for p=1:length(ind_i) % points in image i
        if group_out( ind_i(p) )~=0 % if the point has been classified
            plot(Xi(p,1),Xi(p,2),'.','Color',colors(group_out( ind_i(p) ),:),'MarkerSize',15)
        end
    end
    
end



