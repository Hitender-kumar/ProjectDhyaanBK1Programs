%% FINAL GROUP ANALYSIS (EXPERT vs CONTROL)
clear; clc; close all;

base_path = 'C:\Users\VAIBHAV\Desktop\NSP_GrantProposal\ProjectDhyaanBK1Programs-master\data\ftDataAvgRef\';

%% =========================
% SUBJECT PAIRS + GENDER
%% =========================

pairs = {
    '019CKa','022SSP','M'
    '096MS','026HM','M'
    '040VS','100UK','M'
    '012GK','093AK','M'
    '095KM','075AD','M'
    '056PR','086AB','F'
    '052PR','082MS','F'
    '059MS','102AS','F'
    '013AR','064PK','F'
    '074KS','084AK','F'
};

nPairs = size(pairs,1);

%% =========================
% CHANNELS + BANDS
%% =========================

occipital_labels = {'O1','Oz','O2','POz','PO3','PO4','PO7','PO8','Iz'};
bands = [1 4; 4 8; 8 12; 13 30];

%% =========================
% FUNCTIONS
%% =========================

function feats = extract_features(file_name, occipital_labels, bands)
    data = load(file_name);
    data_struct = data.data;

    fs = data_struct.fsample;
    idx = find(ismember(data_struct.label, occipital_labels));

    feats = [];

    for t = 1:length(data_struct.trial)
        trial = data_struct.trial{t};
        signal = mean(trial(idx,:),1);

        [pxx,f] = pwelch(signal,[],[],[],fs);

        feat = zeros(1,size(bands,1));
        for b = 1:size(bands,1)
            band_idx = f>=bands(b,1) & f<=bands(b,2);
            feat(b) = mean(pxx(band_idx));
        end

        feats = [feats; feat];
    end

    feats = log(feats + 1e-8);
    feats = zscore(feats);
end

function sim = compute_similarity(X1,X2)
    vals = [];
    for i=1:size(X1,1)
        for j=1:size(X2,1)
            vals = [vals; abs(corr(X1(i,:)',X2(j,:)'))];
        end
    end
    sim = mean(vals);
end

function [CII_base, CII_med] = compute_CII(subj_path, occipital_labels, bands)

    EO1 = extract_features(fullfile(subj_path,'EO1_ep_v8.mat'), occipital_labels, bands);
    EC1 = extract_features(fullfile(subj_path,'EC1_ep_v8.mat'), occipital_labels, bands);
    G1  = extract_features(fullfile(subj_path,'G1_ep_v8.mat'), occipital_labels, bands);

    EO2 = extract_features(fullfile(subj_path,'EO2_ep_v8.mat'), occipital_labels, bands);
    EC2 = extract_features(fullfile(subj_path,'EC2_ep_v8.mat'), occipital_labels, bands);
    G2  = extract_features(fullfile(subj_path,'G2_ep_v8.mat'), occipital_labels, bands);

    M1  = extract_features(fullfile(subj_path,'M1_ep_v8.mat'), occipital_labels, bands);
    M2  = extract_features(fullfile(subj_path,'M2_ep_v8.mat'), occipital_labels, bands);

    % Baseline
    b = [ ...
        compute_similarity(EO1,EC1), ...
        compute_similarity(EO1,G1), ...
        compute_similarity(EC1,G1)];

    CII_base = mean(b);

    % Meditation
    m = [ ...
        compute_similarity(M1,EO1), ...
        compute_similarity(M1,EC1), ...
        compute_similarity(M1,G1), ...
        compute_similarity(M2,EO2), ...
        compute_similarity(M2,EC2), ...
        compute_similarity(M2,G2)];

    CII_med = mean(m);
end

%% =========================
% MAIN LOOP
%% =========================

expert_base = []; expert_med = [];
control_base = []; control_med = [];

gender_labels = {};

for i = 1:nPairs

    exp_id = pairs{i,1};
    ctrl_id = pairs{i,2};
    gender = pairs{i,3};

    fprintf('\nPair %d: %s vs %s\n', i, exp_id, ctrl_id);

    [b1,m1] = compute_CII(fullfile(base_path,exp_id), occipital_labels, bands);
    [b2,m2] = compute_CII(fullfile(base_path,ctrl_id), occipital_labels, bands);

    expert_base = [expert_base; b1];
    expert_med  = [expert_med; m1];

    control_base = [control_base; b2];
    control_med  = [control_med; m2];

    gender_labels{i} = gender;
end

%% =========================
% GROUP RESULTS
%% =========================

fprintf('\n===== GROUP MEANS =====\n');

fprintf('Expert Baseline: %.4f\n', mean(expert_base));
fprintf('Expert Meditation: %.4f\n', mean(expert_med));
fprintf('Control Baseline: %.4f\n', mean(control_base));
fprintf('Control Meditation: %.4f\n', mean(control_med));

%% =========================
% STATISTICAL TESTS
%% =========================

% Paired Expert vs Control
[p_base,~] = ttest(expert_base, control_base);
[p_med,~]  = ttest(expert_med, control_med);

% Within-group (state effect)
[p_exp_state,~] = ttest(expert_med, expert_base);
[p_ctrl_state,~] = ttest(control_med, control_base);

fprintf('\n===== P-VALUES =====\n');
fprintf('Baseline Expert vs Control: %.6e\n', p_base);
fprintf('Meditation Expert vs Control: %.6e\n', p_med);
fprintf('Expert Baseline vs Meditation: %.6e\n', p_exp_state);
fprintf('Control Baseline vs Meditation: %.6e\n', p_ctrl_state);

%% =========================
% EFFECT SIZE (Cohen''s d)
%% =========================

d_base = mean(expert_base - control_base) / std(expert_base - control_base);
d_med  = mean(expert_med - control_med) / std(expert_med - control_med);

fprintf('\n===== EFFECT SIZE =====\n');
fprintf('Cohen d (Baseline): %.4f\n', d_base);
fprintf('Cohen d (Meditation): %.4f\n', d_med);

%% =========================
% IMPROVEMENT (INTERACTION)
%% =========================

delta_exp = expert_med - expert_base;
delta_ctrl = control_med - control_base;

[p_interaction,~] = ttest(delta_exp, delta_ctrl);

fprintf('\n===== INTERACTION =====\n');
fprintf('Expert Improvement: %.4f\n', mean(delta_exp));
fprintf('Control Improvement: %.4f\n', mean(delta_ctrl));
fprintf('Interaction p-value: %.6e\n', p_interaction);

%% =========================
% GENDER ANALYSIS
%% =========================

male_idx = strcmp(gender_labels,'M');
female_idx = strcmp(gender_labels,'F');

fprintf('\n===== GENDER ANALYSIS =====\n');

fprintf('Male Expert Mean: %.4f\n', mean(expert_med(male_idx)));
fprintf('Female Expert Mean: %.4f\n', mean(expert_med(female_idx)));

%% =========================
% PLOTS
%% =========================

% Main comparison
figure;
boxplot([expert_base, control_base, expert_med, control_med], ...
    'Labels', {'Exp-Base','Ctrl-Base','Exp-Med','Ctrl-Med'});
ylabel('CII');
title('Expert vs Control Comparison');
grid on;

% Paired difference plot
figure;
plot(expert_base - control_base, '-o'); hold on;
plot(expert_med - control_med, '-o');
legend('Baseline','Meditation');
title('Paired Differences (Expert - Control)');
ylabel('CII Difference');
grid on;

% Improvement comparison
figure;
boxplot([delta_exp, delta_ctrl], 'Labels', {'Expert','Control'});
ylabel('CII Change');
title('Meditation Effect (Improvement)');
grid on;