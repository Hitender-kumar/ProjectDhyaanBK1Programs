%% CII USING PPC (TRIAL-LEVEL SIMILARITY)
clear; clc; close all;

%% =========================
% PATH
%% =========================

base_path = 'D:\ftavgdata\013AR';

files = {
    'EC1_ep_v8.mat'
    'EO1_ep_v8.mat'
    'G1_ep_v8.mat'
};

labels_plot = {'EC','EO','GR'};

%% =========================
% CHANNELS
%% =========================

occipital_labels = {'O1','Oz','O2','POz','PO3','PO4','PO7','PO8','Iz'};
frontal_labels   = {'Fz','F1','F2','F3','F4'};

%% =========================
% BANDS (for PPC)
%% =========================

bands = [
    4 8;    % theta
    8 12;   % alpha
    13 30   % beta
];

nBands = size(bands,1);

%% =========================
% FEATURE EXTRACTION (PPC)
%% =========================

features = cell(3,1);

for c = 1:3
    
    file_name = fullfile(base_path, files{c});
    data = load(file_name);
    data_struct = data.data;
    
    fs = data_struct.fsample;
    
    occ_idx = find(ismember(data_struct.label, occipital_labels));
    fro_idx = find(ismember(data_struct.label, frontal_labels));
    
    feats = [];
    
    for t = 1:length(data_struct.trial)
        
        trial = data_struct.trial{t};
        
        sig_occ = mean(trial(occ_idx,:),1);
        sig_fro = mean(trial(fro_idx,:),1);
        
        % High-pass
        [b,a] = butter(3, 1/(fs/2), 'high');
        sig_occ = filtfilt(b,a,sig_occ);
        sig_fro = filtfilt(b,a,sig_fro);
        
        feat = zeros(1,nBands);
        
        for bnd = 1:nBands
            
            band = bands(bnd,:);
            
            % Bandpass
            s1 = bandpass_filt(sig_fro, fs, band);
            s2 = bandpass_filt(sig_occ, fs, band);
            
            % Phase
            ph1 = angle(hilbert(s1));
            ph2 = angle(hilbert(s2));
            
            dphi = ph1 - ph2;
            
            % PPC
            feat(bnd) = compute_ppc(dphi);
        end
        
        feats = [feats; feat]; %#ok<AGROW>
    end
    
    % Normalize
    feats = zscore(feats);
    
    features{c} = feats;
end

%% =========================
% COMPUTE CROSS-SIMILARITY
%% =========================

sim_matrix = zeros(3);

for i = 1:3
    for j = 1:3
        
        Xi = features{i};
        Xj = features{j};
        
        sims = [];
        
        for a = 1:size(Xi,1)
            for b = 1:size(Xj,1)
                sims = [sims; corr(Xi(a,:)', Xj(b,:)')]; %#ok<AGROW>
            end
        end
        
        sim_matrix(i,j) = mean(sims);
    end
end

%% =========================
% DISPLAY
%% =========================

disp('=== PPC SIMILARITY MATRIX ===');
disp(sim_matrix);

figure;
imagesc(sim_matrix);
colorbar;
title('PPC Trial-Level Similarity Matrix');
xticks(1:3); yticks(1:3);
xticklabels(labels_plot);
yticklabels(labels_plot);

%% =========================
% TRUE CII
%% =========================

CII = mean(sim_matrix(~eye(3)));

fprintf('\nPPC CII: %.3f\n', CII);

%% =========================
% FUNCTIONS
%% =========================

function y = bandpass_filt(x, fs, band)
    [b,a] = butter(3, band/(fs/2), 'bandpass');
    y = filtfilt(b,a,x);
end

function ppc = compute_ppc(dphi)
    dphi = dphi(:);
    N = length(dphi);
    z = exp(1i * dphi);
    ppc = (abs(sum(z))^2 - N) / (N*(N-1));
    ppc = real(ppc);
end