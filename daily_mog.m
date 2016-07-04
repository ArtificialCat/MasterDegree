clear all
clc

%% daily test


load('141026_data.mat');        %% normalized data
%   load('141106_data.mat');        %% original data

Date=unique(M(:,end));
win_size = 3;
mog_auc_mean = zeros(size(Date,1)-win_size,1);
parzen_auc_mean = zeros(size(Date,1)-win_size,1);

for date_idx=1:size(Date,1)-win_size
    
    trn_idx = ismember(M(:,end),Date(date_idx:date_idx+win_size,1));
    tst_idx = ismember(M(:,end),Date(date_idx+win_size+1:date_idx+2*win_size,1));
    
    
    trn_inst = M(trn_idx(:,1)==1,1:end-2);
    trn_label = M(trn_idx(:,1)==1,end-1);
    tst_inst = M(tst_idx(:,1)==1,1:end-2);
    tst_label = M(tst_idx(:,1)==1,end-1);
    
    temp = prdataset(trn_inst,trn_label);
    temp2 = prdataset(tst_inst,tst_label);

    trn = oc_set(temp,'0');
    tst = oc_set(temp2,'0');

    w = mog_dd(trn,0.01);

    mog_auc = dd_auc(tst*w*dd_roc);
    mog_auc_mean(date_idx,1) = mog_auc;
 
    w2= parzen_dd(trn,0.1);
    parzen_auc = dd_auc(tst*w2*dd_roc);
    parzen_auc_mean(date_idx,1) = parzen_auc;
end


%% raw data test

% clear;
% clc;
% 
% load '20141103_data.mat'            %% raw data(normalize: x, feature selection: x)
% 
% NumBag=5;           %% fold ¼ö
% theta = 0.05;       %% Theta °ª (Mixture of gaussian)
% 
% temp = prdataset(M(:,1:end-1),M(:,end));
% data = oc_set(temp,'0');
% auc = zeros(1,NumBag);
% I=NumBag;
% 
% 
% for i=1:NumBag
%     
%     [trn,tst,I] = dd_crossval(data,I);
%     w = mog_dd(trn,theta);
%     auc(1,i) = dd_auc(tst*w*dd_roc);
% end
% 
% auc_mean = mean(auc);






