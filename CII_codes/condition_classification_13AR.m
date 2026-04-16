%% CONDITION CLASSIFICATION - SUBJECT 13AR
clear; clc; close all;

%% =========================
% PATH SETTINGS (EDIT THIS)
%% =========================

base_path = 'C:\Users\VAIBHAV\Desktop\NSP_GrantProposal\ProjectDhyaanBK1Programs-master\data\ftDataAvgRef\013AR\';

files = {
    'EC1_ep_v8.mat'   % 1
    'EO1_ep_v8.mat'   % 2
    'G1_ep_v8.mat'    % 3
};

labels_plot = {'EC','EO','GR'};

%% =========================
% CHANNELS + BANDS
%% =========================

occipital_labels = {'O1','Oz','O2','POz','PO3','PO4','PO7','PO8','Iz'};

bands = [
    1 4;    % delta
    4 8;    % theta
    8 12;   % alpha
    13 30   % beta
];

nBands = size(bands,1);

%% =========================
% FEATURE EXTRACTION
%% =========================

X = [];
y = [];

for c = 1:length(files)
    
    file_name = fullfile(base_path, files{c});
    
    data = load(file_name);
    data_struct = data.data;
    
    fs = data_struct.fsample;
    nTrials = length(data_struct.trial);
    
    occipital_idx = find(ismember(data_struct.label, occipital_labels));
    
    fprintf('Processing %s | Trials: %d\n', files{c}, nTrials);
    
    for t = 1:nTrials
        
        trial = data_struct.trial{t};
        signal = mean(trial(occipital_idx,:),1);
        
        [pxx, f] = pwelch(signal, [], [], [], fs);
        
        feat = zeros(1, nBands);
        
        for b = 1:nBands
            idx = f >= bands(b,1) & f <= bands(b,2);
            feat(b) = mean(pxx(idx));
        end
        
        X = [X; feat]; %#ok<AGROW>
        y = [y; c]; %#ok<AGROW>
    end
end

%% =========================
% LOG + NORMALIZATION
%% =========================

X = log(X + 1e-8);
X = zscore(X);

%% =========================
% CLASSIFICATION (5-FOLD CV)
%% =========================

cv = cvpartition(y, 'KFold', 5);

acc = zeros(cv.NumTestSets,1);

for k = 1:cv.NumTestSets
    
    train_idx = cv.training(k);
    test_idx  = cv.test(k);
    
    mdl = fitcecoc(X(train_idx,:), y(train_idx));
    
    y_pred = predict(mdl, X(test_idx,:));
    
    acc(k) = mean(y_pred == y(test_idx));
end

accuracy = mean(acc);

fprintf('\nClassification Accuracy: %.3f\n', accuracy);

%% =========================
% CONFUSION MATRIX
%% =========================

mdl = fitcecoc(X, y);
y_pred_all = predict(mdl, X);

figure;
confusionchart(y, y_pred_all);
title('Confusion Matrix (All Data)');