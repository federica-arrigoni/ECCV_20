
function [missrate_known,known_rate,group_out,ind_known]=compute_missrate(group_out,labels_gt)

% Classified points
ind_known=find(group_out~=0);

% Alignment (only over classified points)
group_out(ind_known)=bestMap(labels_gt(ind_known),group_out(ind_known));

% compute misclassification rate
missrate_known = sum(labels_gt(ind_known) ~= group_out(ind_known)) / length(ind_known);
known_rate=length(ind_known)/length(labels_gt);


end
