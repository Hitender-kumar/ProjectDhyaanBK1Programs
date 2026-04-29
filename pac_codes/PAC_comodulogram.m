%% ================================================================
%  PAC Comodulogram single subject single condtition

clear; clc;

%% ===== USER SETTINGS =====

subj_folder = 'D:\ftavgdata\013AR';
cond_file   = 'M2_ep_v8.mat';

n_surr = 50;   

phase_centers = 5:1:18;
phase_bw      = 2;

amp_centers   = 10:5:100;
amp_bw        = 10;

%% ===== FIXED PARAMETERS =====

fs     = 1000;
n_bins = 18;

% Occipital (1-indexed): O1=16, Oz=17, O2=18
occ_idx   = [16, 17, 18];

% Frontal (1-indexed): Fz=2, F3=3, F4=29, AFz=35, AF3=34, AF4=62, F1=36, F2=63
front_idx = [2, 3, 29, 35, 34, 62, 36, 63];

%% ===== LOAD DATA =====

fprintf('Loading: %s\n', fullfile(subj_folder, cond_file));
raw = load(fullfile(subj_folder, cond_file));
ft  = raw.data;

trials_cell = ft.trial;
nTrials     = length(trials_cell);
bad_elecs   = ft.badElecs(:);

fprintf('Trials: %d | Bad elecs: %s\n', nTrials, mat2str(bad_elecs'));

%% ===== CHANNEL SELECTION =====

occ_clean   = setdiff(occ_idx,   bad_elecs);
front_clean = setdiff(front_idx, bad_elecs);

fprintf('Occipital channels used : %d\n', numel(occ_clean));
fprintf('Frontal channels used   : %d\n\n', numel(front_clean));

if isempty(occ_clean),   error('All occipital channels are bad!'); end
if isempty(front_clean), error('All frontal channels are bad!');   end

%% ===== BUILD PER-TRIAL SIGNALS =====

[b_hp, a_hp] = butter(3, 1/(fs/2), 'high');

occ_trials   = cell(1, nTrials);
front_trials = cell(1, nTrials);

for t = 1:nTrials
    occ_sig   = mean(ft.trial{t}(occ_clean,   :), 1);
    front_sig = mean(ft.trial{t}(front_clean, :), 1);
    occ_trials{t}   = filtfilt(b_hp, a_hp, occ_sig);
    front_trials{t} = filtfilt(b_hp, a_hp, front_sig);
end

%% ===== COMPUTE BOTH COMODULOGRAM DIRECTIONS =====

nP = numel(phase_centers);
nA = numel(amp_centers);
 

MI_frt_occ  = zeros(nP, nA);
Zsc_frt_occ = zeros(nP, nA);
MI_occ_frt  = zeros(nP, nA);
Zsc_occ_frt = zeros(nP, nA);

fprintf('=== Cross-Regional Comodulogram ===\n');
fprintf('    Direction A: Frontal phase -> Occ amplitude (top-down)\n');
fprintf('    Direction B: Occ phase -> Frontal amplitude (bottom-up)\n');
fprintf('    %d x %d = %d pairs, %d surrogates/trial\n\n', ...
        nP, nA, nP*nA, n_surr);

tic;

for pi = 1:nP

    pf    = phase_centers(pi);
    plo   = pf - phase_bw;
    phi_f = pf + phase_bw;

    for ai = 1:nA

        af  = amp_centers(ai);
        alo = max(phi_f + 2, af - amp_bw);
        ahi = min(fs/2 - 1, af + amp_bw);

        if alo >= ahi, continue; end

        % Pre-allocate surrogate arrays
        surr_A = zeros(1, nTrials * n_surr);
        surr_B = zeros(1, nTrials * n_surr);
        cnt    = 0;

        trial_MI_A = zeros(1, nTrials);
        trial_MI_B = zeros(1, nTrials);

        for t = 1:nTrials

            n_t     = numel(occ_trials{t});

            
            ph_frt  = angle(hilbert(do_bandpass(front_trials{t}, plo, phi_f, fs)));
            ph_occ  = angle(hilbert(do_bandpass(occ_trials{t},   plo, phi_f, fs)));

            
            env_occ = abs(hilbert(do_bandpass(occ_trials{t},   alo, ahi, fs)));
            env_frt = abs(hilbert(do_bandpass(front_trials{t}, alo, ahi, fs)));

          
            trial_MI_A(t) = tort_MI(ph_frt, env_occ, n_bins);   
            trial_MI_B(t) = tort_MI(ph_occ, env_frt, n_bins);   

            for sv = 1:n_surr
                shift = randi([round(n_t*0.1), round(n_t*0.9)]);
                cnt   = cnt + 1;
                surr_A(cnt) = tort_MI(ph_frt, circshift(env_occ, shift), n_bins);
                surr_B(cnt) = tort_MI(ph_occ, circshift(env_frt, shift), n_bins);
            end

        end 

        real_A = mean(trial_MI_A);
        real_B = mean(trial_MI_B);
        s_used = surr_A(1:cnt);

        MI_frt_occ(pi,ai)  = real_A;
        Zsc_frt_occ(pi,ai) = (real_A - mean(surr_A(1:cnt))) / (std(surr_A(1:cnt)) + eps);

        MI_occ_frt(pi,ai)  = real_B;
        Zsc_occ_frt(pi,ai) = (real_B - mean(surr_B(1:cnt))) / (std(surr_B(1:cnt)) + eps);

    end 
    fprintf('  Phase %2d-%2d Hz | %d/%d | %.1f min\n', ...
            plo, phi_f, pi, nP, toc/60);

end
fprintf('\nDone in %.1f minutes.\n\n', toc/60);

%% ===== FIGURE 1: Frt->Occ Z-score (top-down) =====

max_zA = max(Zsc_frt_occ(:));

figure('Name', 'Comodulogram: Frontal->Occipital (Top-down)', ...
       'Position', [50 50 750 520]);
pcolor(amp_centers, phase_centers, Zsc_frt_occ);
shading flat; colormap(hot);
clim([0, max(max_zA * 1.1, 3)]);
cb = colorbar; cb.Label.String = sprintf('Z-score (n=%d surr/trial)', n_surr);
xlabel('Amplitude Frequency (Hz)', 'FontSize', 12);
ylabel('Phase Frequency (Hz)',     'FontSize', 12);
title(sprintf('TOP-DOWN: Frontal phase \\rightarrow Occ amplitude\n%s | %s | max Z=%.2f', ...
      strrep(cond_file,'_ep_v8.mat',''), subj_folder, max_zA), ...
      'FontSize', 10, 'Interpreter', 'tex');
hold on;
[~,hc] = contour(amp_centers, phase_centers, Zsc_frt_occ,[1.96 1.96],'c-');
set(hc,'LineWidth',2,'DisplayName','Z=1.96');
patch([30 80 80 30 30],[8 8 12 12 8],'w','FaceAlpha',0,'EdgeColor','w','LineWidth',2);
text(55,10,'\alpha_{Frt}\rightarrow\gamma_{Occ}','Color','w','FontWeight','bold', ...
     'FontSize',10,'HorizontalAlignment','center');
legend(hc,'Z=1.96 (p<0.05)','Location','northeast');
hold off;

%% ===== FIGURE 2: Occ->Frt Z-score (bottom-up) =====

max_zB = max(Zsc_occ_frt(:));
    
figure('Name', 'Comodulogram: Occipital->Frontal (Bottom-up)', ...
       'Position', [100 100 750 520]);
pcolor(amp_centers, phase_centers, Zsc_occ_frt);
shading flat; colormap(hot);
clim([0, max(max_zB * 1.1, 3)]);
cb = colorbar; cb.Label.String = sprintf('Z-score (n=%d surr/trial)', n_surr);
xlabel('Amplitude Frequency (Hz)', 'FontSize', 12);
ylabel('Phase Frequency (Hz)',     'FontSize', 12);
title(sprintf('BOTTOM-UP: Occ phase \\rightarrow Frontal amplitude\n%s | %s | max Z=%.2f', ...
      strrep(cond_file,'_ep_v8.mat',''), subj_folder, max_zB), ...
      'FontSize', 10, 'Interpreter', 'tex');
hold on;
[~,hc] = contour(amp_centers, phase_centers, Zsc_occ_frt,[1.96 1.96],'c-');
set(hc,'LineWidth',2,'DisplayName','Z=1.96');
patch([30 80 80 30 30],[8 8 12 12 8],'w','FaceAlpha',0,'EdgeColor','w','LineWidth',2);
text(55,10,'\alpha_{Occ}\rightarrow\gamma_{Frt}','Color','w','FontWeight','bold', ...
     'FontSize',10,'HorizontalAlignment','center');
legend(hc,'Z=1.96 (p<0.05)','Location','northeast');
hold off;



%% ===== LOCAL FUNCTIONS =====

function sig_filt = do_bandpass(sig, lo, hi, fs)
    lo = max(lo, 2);
    hi = min(hi, fs/2 - 1);
    if lo >= hi
        sig_filt = zeros(size(sig));
        return;
    end
    [b, a]   = butter(3, [lo hi]/(fs/2), 'bandpass');
    sig_filt = filtfilt(b, a, sig);
end

function MI = tort_MI(phase, amp, n_bins)
    edges      = linspace(-pi, pi, n_bins + 1);
    amp_binned = zeros(1, n_bins);
    for k = 1:n_bins
        idx = phase >= edges(k) & phase < edges(k+1);
        if any(idx)
            amp_binned(k) = mean(amp(idx));
        end
    end
    s = sum(amp_binned) + eps;
    amp_dist = amp_binned / s;
    H  = -sum(amp_dist .* log(amp_dist + eps));
    MI = (log(n_bins) - H) / log(n_bins);
end