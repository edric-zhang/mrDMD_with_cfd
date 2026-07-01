close all;
clc;
clear;

refresh = false;
if refresh || ~exist('X', 'var')
    fprintf('Loading cylinder_2D_uv_ONE_MIDDLE_Z_cropped_downsampled.mat...\n');
    load('cylinder_2D_uv_ONE_MIDDLE_Z_cropped_downsampled.mat');
end

X = single(X);

if ~exist('npoints', 'var')
    if mod(size(X,1), 2) ~= 0
        error('Expected X to be stacked as [u; v], but size(X,1) is odd.');
    end
    npoints = size(X,1) / 2;
end

if size(X,1) ~= 2*npoints
    error('Expected size(X,1) to equal 2*npoints for stacked [u; v] snapshots.');
end

if ~exist('x', 'var') || ~exist('y', 'var')
    error('The CFD dataset should include point coordinates x and y.');
end

x = x(:);
y = y(:);
if numel(x) ~= npoints || numel(y) ~= npoints
    error('Expected x and y to each contain npoints entries.');
end

if exist('dt', 'var')
    cfd_dt = dt;
else
    cfd_dt = 1; % Use snapshot index as time if the MAT file does not include physical dt.
end

X_all = X;
if size(X_all, 2) >= 401
    X_test = single(X_all(:,401:end));
else
    X_test = single([]);
end

% Keep the same 1:400 working window used by the original sweep.
% Increase this later once you want the sweep to exploit more CFD snapshots.
X = X_all(:, 1:min(400, size(X_all, 2)));
clear X_all;

%% Setting the parameters
X_total_size = size(X,2);
frame_start = 1;
frame_end   = X_total_size;
X = X(:, frame_start:frame_end);
[n, m] = size(X);
numsnapshots = size(X, 2);   % keep this as full working snapshot length

dt = cfd_dt;
total_time = dt * (m - 1);               % True total duration of analyzed snapshots
t = (0:dt:total_time)';
%% STARTING MRDMD
L = 3; 
maxJ = 2^(L-1); 
matrices = cell(1, maxJ);       % List of matrices that we update to run DMD through
matrices{1} = X;

% Tracking variables
list_ml = zeros(L, maxJ);           % List of truncation number for each level, bin
list_b = cell(L, maxJ);             % List of starting vector coefficients for each level, bin mode
list_w = cell(L, maxJ);             % List of eigenvalues for each level, bin mode
list_modes = cell(L, maxJ);         % List of modes for each level, bin
list_bin_widths = zeros(L, maxJ);   % List of bin widths for each level, bin
list_t_start = zeros(L, maxJ);      % List of bin start points for each level, bin
list_t_start(1,1) = 1;
list_bin_widths(1,1) = m;
freq_threshold_hz = 1000;

for i = 1:L
    J = 2^(i-1);
    next_matrices = cell(1, 2*J);   % Setting an empty list of matrices for the NEXT level
    count = 1;                      % Used to go to the next spot in the next_matrices list
    
    for j = 1:J
        A = matrices{j};            % Get the next matrix from the previous matrices list
        
        if isempty(A) || size(A, 2) < 10 || any(isnan(A(:)))    % If A is empty or too small, 
            next_matrices{count} = []; count = count + 1;       % then add an empty matrix for the 
            next_matrices{count} = []; count = count + 1;       % next list and skip the process
            continue;
        end
        
        current_start = list_t_start(i, j);                     % For (1,1), starts at 0
        current_width = size(A, 2);
        list_bin_widths(i, j) = current_width;                  % For (1,1), starts at 0
        
        level_dt = dt;                                          % CFD dt, or snapshot-index time if dt was not in the MAT file
        [modes, D, b] = dmd(A);                                 % Run DMD
        timeperiod = size(D, 1);                                % Size of D (depends on our r value)
        
        
        % Calculate continuous frequencies
        freqs = zeros(timeperiod, 1);                   
        for modenum = 1:timeperiod
            lambda = D(modenum, modenum);                       % For each number in our number of modes
            omega = log(lambda) / level_dt;                     % Making this a continuous eigenvalue
            freqs(modenum) = abs(imag(omega)) / (2 * pi);       % Finding the corresponding frequency
        end
        

        [sorted_freqs, sort_idx] = sort(freqs);                 % Sort frequencies ascending
        D = D(sort_idx, sort_idx);                              % Sort D, modes, and b in the same way
        modes = modes(:, sort_idx);
        b = b(sort_idx);
        
        ml = find(sorted_freqs >= freq_threshold_hz, 1, 'first'); % Find the first frequency to be higher than the threshold
        
        if isempty(ml)                                          % If none are higher, it flags ALL MODES
            ml = timeperiod + 1;                                % as being SLOW MODES
            slow_inds = 1:timeperiod;
        else
            slow_inds = 1:(ml-1);                               % Or else it flags all up to that index
        end
        
        if ~isempty(slow_inds)                                  % As long as there exist slow modes
            eigs_slow = diag(D(slow_inds, slow_inds));
            mag = abs(eigs_slow);                               
            over = mag > 1.0;                                   % If the magnitude of the discrete eigenvalues
            eigs_slow(over) = eigs_slow(over) ./ mag(over);     % is greater than 1, project it onto 1 unit circle
                                                                
            time_powers = eigs_slow .^ (0:current_width-1);     % basically Omega Matrix
            slowmatrix = modes(:, slow_inds) * (b(slow_inds) .* time_powers);
        else
            slowmatrix = zeros(n, current_width, 'like', X);              % If there are no slow modes, just make zeros
        end
        
        fastmatrix = A - slowmatrix; 
        clear A slowmatrix;
        % Storing all the variables
        list_ml(i, j) = ml;                                     
        list_w{i, j} = diag(D(slow_inds, slow_inds));           % Store eigenvalues as a vector
        list_b{i, j} = b(slow_inds);
        list_modes{i, j} = modes(:, slow_inds);
        
        % Splitting Step: Split the new fastmatrix into half
        midpoint = floor(current_width / 2);
        A1 = fastmatrix(:, 1:midpoint);
        A2 = fastmatrix(:, (midpoint + 1):end);
        clear fastmatrix;
        left_child_idx = 2*j - 1;
        next_matrices{left_child_idx} = A1;
        right_child_idx = 2*j;
        next_matrices{right_child_idx} = A2;
        if i < L
            list_t_start(i+1, left_child_idx) = current_start;             % Setting starts for the NEXT LEVEL
            list_bin_widths(i+1, left_child_idx) = midpoint;

            list_t_start(i+1, right_child_idx) = current_start + midpoint;
            list_bin_widths(i+1, right_child_idx) = current_width - midpoint;
        end
        
    clear A1 A2 modes D b freqs sorted_freqs sort_idx time_powers eigs_slow mag over;
    end
    matrices = next_matrices; % Update the new matrix list
    clear next_matrices;
end
clear matrices next_matrices
clear A1 A2 A fastmatrix slowmatrix;

%% Reconstruction Loop

X_rec = zeros(n, m, 'like', X);
for i = 1:L
    J = 2^(i-1);
    for j = 1:J
        if isempty(list_modes{i, j})
            continue;
        end
        modes = list_modes{i, j};
        eigs_slow = list_w{i, j};
        b = list_b{i, j};
        
        t_start = list_t_start(i, j);                                      % Setting start/end - modes only exist locally
        bin_width = list_bin_widths(i, j);
        
        if bin_width == 0
            continue;
        end
        
        t_end = t_start + bin_width - 1;
        mag = abs(eigs_slow);
        over = mag > 1.0;                                                  % Project eigs>1 to 1
        eigs_slow(over) = eigs_slow(over) ./ mag(over);
        
        time_powers = eigs_slow .^ (0:bin_width-1);
        local_rec = modes * (b .* time_powers);
        
        % Protect against matrix index rounding overflows at tree boundaries
        if t_end > m                                                       % Make overflow over m into m
            t_end = m;
            local_rec = local_rec(:, 1:(t_end - t_start + 1));
        end

        X_rec(:, t_start:t_end) = X_rec(:, t_start:t_end) + local_rec;     % Add on the local mode data
        clear local_rec time_powers modes eigs_slow b mag over;
    end
end
X_rec = real(X_rec); 

%{
%% Quick display of what we got: 
total_modes = 0;
for i = 1:L
    J = 2^(i-1);
    for j = 1:J
        if ~isempty(list_modes{i,j})
            % The number of columns equals the number of modes in this bin
            num_modes_in_bin = size(list_modes{i,j}, 2); 
            total_modes = total_modes + num_modes_in_bin;
            fprintf('Level %d, Bin %d has %d slow modes\n', i, j, num_modes_in_bin);
        end
    end
end
fprintf('--- Total modes across all time windows: %d ---\n', total_modes);
%}

%% Error Graph

mrdmd_error = NaN(1, numsnapshots); 
for k = 1:m
    true_snapshot = X(:, k);
    rec_snapshot  = X_rec(:, k);
    
    % Shift the tracking index to map precisely onto the absolute X_total timeline
    absolute_idx = frame_start + k - 1;
    mrdmd_error(absolute_idx) = norm(true_snapshot - rec_snapshot) / (norm(true_snapshot) + eps);
end

figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.5]);

% Plotting against the absolute indices preserves the alignment layout
plot(mrdmd_error * 100, 'r-', 'LineWidth', 1.5, 'DisplayName', 'MRDMD');
hold on;
xlim([1, numsnapshots]); % Lock the visual frame boundary to the total absolute snapshots
ylim([0 200]);
grid on;
legend('Location', 'best');
xlabel('Snapshot (Absolute Timeline)'); 
ylabel('Error %');
title('Snapshot-wise Reconstruction Performance');
clear X_rec true_snapshot rec_snapshot local_rec time_powers;



%% MULTI-TARGET EXTRACTION WITH ACTIVE LOOKUP MENU
fit_start_idx = 263;
fit_end_idx   = 310;
test_start_idx = 311;
test_end_idx   = 350;
plot_start_idx = 100; 
plot_end_idx   = m;
num_plot_snaps = X_total_size;



% Format: [Level, Bin (j), Mode_Index] (Use 0 as a wildcard for ALL)
% Use 1,0,0 - 2,0,0 - 3,0,0 for full period-specific reconstruction
target_coordinates = [
    %1, 0, 0;   
    %2, 0, 0; 
    %3, 0, 0;
    3, 0, 0;
];

fprintf('\n===========================================================');
fprintf('\n   AVAILABLE COORDINATES FOR FRAMES %d TO %d', plot_start_idx, plot_end_idx);
fprintf('\n===========================================================\n');


available_modes = cell(L,maxJ);

for i = 1:L
    J = 2^(i-1);
    for j = 1:J
        if isempty(list_modes{i, j}), continue; end
        
        t_start   = list_t_start(i, j);
        bin_width = list_bin_widths(i, j);
        t_end     = t_start + bin_width - 1;
        
        % Check if this bin overlaps with our visual playback window
        if t_start <= plot_end_idx && t_end >= plot_start_idx
            num_modes_available = size(list_modes{i, j}, 2);
            available_modes{i,j} = list_modes{i,j};
            fprintf('● Level %d, Bin %d (Frames %d to %d)\n', i, j, t_start, t_end);
            fprintf('  └─ Available Mode Indices: [');
            for m_idx = 1:num_modes_available
                if m_idx == num_modes_available
                    fprintf('%d', m_idx);
                else
                    fprintf('%d, ', m_idx);
                end
            end
            fprintf(']\n\n');
        end
    end
end
fprintf('===========================================================\n\n');




%% WEAK SINDY ON MRDMD TEMPORAL COEFFICIENTS
% Goal:
% Learn d/dt of Level 3 modal amplitudes using Level 1-3 modal amplitudes.
%
% This replaces the old spatial regression:
%   L3 spatial modes = Theta(L1-L2 spatial modes)*Xi
%
% with actual dynamics:
%   d(a_L3)/dt = Theta(a_L1,a_L2,a_L3)*Xi


%% 4. Run WSINDy structure/parameter sweep

fprintf('\nRunning weak SINDy CFD structure sweep on mrDMD modal amplitudes...\n');
tic;

input_levels = [1 2];
target_level = 3;

% Columns:
% lambda1, lambda2, gamma,
% top_input_modes_per_level, top_target_modes,
% max_terms_per_equation, max_quadratic_base_terms
%
% CFD note:
% Your previous full-field error was already low, but L3 percent error was huge.
% So this sweep varies structure more aggressively than lambda alone.
param_grid = [
    % Tight gamma around best structure: 5 inputs, 9 targets, 5 terms, quad 4
    0.004 0.002 0.013   5 9   5 4
    0.006 0.004 0.013   5 9   5 4
    0.008 0.006 0.013   5 9   5 4

    0.004 0.002 0.014   5 9   5 4
    0.006 0.004 0.014   5 9   5 4
    0.008 0.006 0.014   5 9   5 4

    0.004 0.002 0.015   5 9   5 4
    0.006 0.004 0.015   5 9   5 4
    0.008 0.006 0.015   5 9   5 4

    0.004 0.002 0.016   5 9   5 4
    0.006 0.004 0.016   5 9   5 4
    0.008 0.006 0.016   5 9   5 4

    0.004 0.002 0.017   5 9   5 4
    0.006 0.004 0.017   5 9   5 4
    0.008 0.006 0.017   5 9   5 4

    % Check whether 10 targets with fewer terms beats 9 targets
    0.004 0.002 0.013   5 10   5 4
    0.006 0.004 0.013   5 10   5 4
    0.008 0.006 0.013   5 10   5 4

    0.004 0.002 0.014   5 10   5 4
    0.006 0.004 0.014   5 10   5 4
    0.008 0.006 0.014   5 10   5 4

    0.004 0.002 0.015   5 10   5 4
    0.006 0.004 0.015   5 10   5 4
    0.008 0.006 0.015   5 10   5 4

    0.004 0.002 0.016   5 10   5 4
    0.006 0.004 0.016   5 10   5 4
    0.008 0.006 0.016   5 10   5 4

    % Check if even fewer terms improves stability
    0.006 0.004 0.013   5 9   4 4
    0.006 0.004 0.014   5 9   4 4
    0.006 0.004 0.015   5 9   4 4
    0.006 0.004 0.016   5 9   4 4
    0.006 0.004 0.017   5 9   4 4

    0.006 0.004 0.013   5 10   4 4
    0.006 0.004 0.014   5 10   4 4
    0.006 0.004 0.015   5 10   4 4
    0.006 0.004 0.016   5 10   4 4
    0.006 0.004 0.017   5 10   4 4

    % Sanity check: 8 targets with best complexity
    0.006 0.004 0.013   5 8   5 4
    0.006 0.004 0.014   5 8   5 4
    0.006 0.004 0.015   5 8   5 4
    0.006 0.004 0.016   5 8   5 4
    0.006 0.004 0.017   5 8   5 4
];
sweep_results = zeros(size(param_grid,1), 11);

for rr = 1:size(param_grid,1)
    lambda1 = param_grid(rr,1);
    lambda2 = param_grid(rr,2);
    gamma = param_grid(rr,3);
    top_input_modes_per_level = param_grid(rr,4);
    top_target_modes = param_grid(rr,5);
    max_terms_per_equation = param_grid(rr,6);
    max_quadratic_base_terms = param_grid(rr,7);

    fprintf('\n=== RUN %d/%d | lambda1 %.4g | lambda2 %.4g | gamma %.4g | inputs %d | targets %d | maxTerms %d | quadBase %d ===\n', ...
        rr, size(param_grid,1), lambda1, lambda2, gamma, ...
        top_input_modes_per_level, top_target_modes, max_terms_per_equation, max_quadratic_base_terms);

    [mean_l3_err, max_l3_err, mean_full_err, max_full_err] = evaluate_wsindy_parameter_set( ...
        list_w, list_b, list_modes, list_t_start, list_bin_widths, ...
        X, dt, m, n, input_levels, target_level, ...
        fit_start_idx, fit_end_idx, test_start_idx, test_end_idx, ...
        top_input_modes_per_level, top_target_modes, ...
        lambda1, lambda2, gamma, max_terms_per_equation, max_quadratic_base_terms);

    sweep_results(rr,:) = [lambda1, lambda2, gamma, ...
        top_input_modes_per_level, top_target_modes, max_terms_per_equation, max_quadratic_base_terms, ...
        mean_l3_err, max_l3_err, mean_full_err, max_full_err];

    fprintf('RESULT | lambda1 %.4g | lambda2 %.4g | gamma %.4g | inputs %d | targets %d | maxTerms %d | quadBase %d | L3 mean %.2f%% | L3 max %.2f%% | full mean %.2f%% | full max %.2f%%\n', ...
        lambda1, lambda2, gamma, top_input_modes_per_level, top_target_modes, ...
        max_terms_per_equation, max_quadratic_base_terms, ...
        mean_l3_err, max_l3_err, mean_full_err, max_full_err);
end

fprintf('\n================ CFD STRUCTURE SWEEP RESULTS ================\n');
result_table = array2table(sweep_results, ...
    'VariableNames', {'lambda1','lambda2','gamma','top_inputs','top_targets','max_terms','quad_base','mean_L3_error','max_L3_error','mean_full_error','max_full_error'});
disp(result_table);

toc;

%% FUNCTIONS
function [modes, D, b] = dmd(X)
    X1 = X(:, 1:end-1);
    Y = X(:, 2:end);
    [U, S, V] = svds(X1, 25);
    %{
    thresh = 1e-8 * sing_vals(1);
    r = sum(sing_vals > thresh);
    r = min([25, r, size(U, 2)]); 
    %}
    r = 25;
    if r == 0, r = 1; end
    U = U(:, 1:r);
    S = S(1:r, 1:r);
    V = V(:, 1:r);
    
    A = (U' * Y * V) / S;
    [W, D] = eig(A);
    
    modes = single(Y * V * (S \ W));
    b = single(pinv(modes) * X1(:, 1)); 
end

function [w_second, labels_second, mode_labels, target_cols, mu_xobs, sigma_xobs] = run_mrdmd_wsindy( ...
    list_w, list_b, list_modes, list_t_start, list_bin_widths, ...
    dt, m, input_levels, target_level, plot_start_idx, plot_end_idx, ...
    top_input_modes_per_level, top_target_modes, lambda1, lambda2, gamma, ...
    max_terms_per_equation, max_quadratic_base_terms)

    if nargin < 14 || isempty(lambda1), lambda1 = 0.006; end
    if nargin < 15 || isempty(lambda2), lambda2 = 0.004; end
    if nargin < 16 || isempty(gamma), gamma = 1e-2; end
    if nargin < 17 || isempty(max_terms_per_equation), max_terms_per_equation = 8; end
    if nargin < 18 || isempty(max_quadratic_base_terms), max_quadratic_base_terms = 5; end

    m_interval = plot_end_idx - plot_start_idx + 1;

    xobs = [];
    mode_labels = {};
    mode_levels = [];
    Lmax = size(list_modes, 1);

    for lev = 1:Lmax
        J = 2^(lev-1);

        for bin = 1:J
            if isempty(list_modes{lev, bin})
                continue;
            end

            eigs_slow = list_w{lev, bin};
            b = list_b{lev, bin};
            bin_width = list_bin_widths(lev, bin);
            t_start = list_t_start(lev, bin);
            t_end = t_start + bin_width - 1;

            if bin_width <= 1 || t_end > m
                continue;
            end

            if t_start > plot_end_idx || t_end < plot_start_idx
                continue;
            end

            local_start = max(t_start, plot_start_idx);
            local_end = min(t_end, plot_end_idx);

            rel_start = local_start - t_start;
            rel_end = local_end - t_start;

            insert_start = local_start - plot_start_idx + 1;
            insert_end = local_end - plot_start_idx + 1;

            for mode_idx = 1:length(eigs_slow)
                a_local = b(mode_idx) * eigs_slow(mode_idx).^(rel_start:rel_end);

                a_full = zeros(m_interval, 1);
                a_full(insert_start:insert_end) = real(a_local(:));

                xobs = [xobs, a_full];
                mode_labels{end+1} = sprintf('L%d B%d Mode%d', lev, bin, mode_idx);
                mode_levels(end+1) = lev;
            end
        end
    end

    tobs = (0:m_interval-1)' * dt;

    keep_cols = false(1, length(mode_levels));

    for lev = input_levels
        lev_cols = find(mode_levels == lev);
        if isempty(lev_cols)
            continue;
        end

        lev_energy = rms(xobs(:, lev_cols), 1);
        [~, ord] = sort(lev_energy, 'descend');
        keep_count = min(top_input_modes_per_level, length(lev_cols));
        keep_cols(lev_cols(ord(1:keep_count))) = true;
    end

    target_level_cols = find(mode_levels == target_level);
    target_energy = rms(xobs(:, target_level_cols), 1);
    [~, target_ord] = sort(target_energy, 'descend');
    keep_target_count = min(top_target_modes, length(target_level_cols));
    keep_cols(target_level_cols(target_ord(1:keep_target_count))) = true;

    xobs = xobs(:, keep_cols);
    mode_labels = mode_labels(keep_cols);
    mode_levels = mode_levels(keep_cols);
    target_cols = find(mode_levels == target_level);

    fprintf('\nReduced WSINDy library to %d modes: %d input modes, %d target modes.\n', ...
        size(xobs, 2), sum(ismember(mode_levels, input_levels)), length(target_cols));

    weights = [];
    polys = 1;
    trigs = [];
    custom_tags = [];
    custom_fcns = {};
    phi_class = 1;
    max_d = 1;
    tau = 10^-16;
    tauhat = -1;
    K_frac = min(50, floor(m_interval/2));
    overlap_frac = 1;
    relax_AG = 0;
    scale_Theta = 2;
    useGLS = 0;

    include_l3_peer_coupling = true;
    include_quadratic_terms = true;
    include_cross_terms = true;
    max_active_inputs_per_equation = 4;
    max_peer_l3_per_equation = 2;
    relative_term_tol = 0.03;
    max_nonlinear_terms_per_equation = 4;

    alpha_loss = 0.8;
    overlap_frac_ag = 0.8;
    mt_ag_fac = [1;0];
    pt_ag_fac = [1;0];
    run_sindy = 0;
    useFD = 2;
    smoothing_window = 0;

    mu_xobs = mean(xobs, 1);
    sigma_xobs = std(xobs, 0, 1);
    sigma_xobs(sigma_xobs == 0) = 1;

    xobs = (xobs - mu_xobs) ./ sigma_xobs;

    [w_linear,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
        wsindy_ode_fun_capped(xobs,tobs,weights, ...
        polys,trigs,custom_tags,custom_fcns, ...
        phi_class,max_d,tau,tauhat,K_frac,overlap_frac,relax_AG, ...
        scale_Theta,useGLS,lambda1,gamma,alpha_loss, ...
        overlap_frac_ag,pt_ag_fac,mt_ag_fac,run_sindy,useFD,smoothing_window);

    w_second = cell(length(target_cols),1);
    labels_second = cell(length(target_cols),1);

    n_base = size(xobs,2);

    for kk = 1:length(target_cols)

        target_col = target_cols(kk);

        coef1 = w_linear(:, target_col);
        tol1 = 1e-6 * max(abs(coef1));

        active_linear = find(abs(coef1(1:n_base)) > tol1);

        active_inputs = active_linear(ismember(mode_levels(active_linear), input_levels));

        if ~isempty(active_inputs)
            [~, ord] = sort(abs(coef1(active_inputs)), 'descend');
            active_inputs = active_inputs(ord(1:min(max_active_inputs_per_equation, length(active_inputs))));
        end

        active_peers = [];
        if include_l3_peer_coupling
            active_peers = active_linear(ismember(mode_levels(active_linear), target_level));
            active_peers(active_peers == target_col) = [];

            if ~isempty(active_peers)
                [~, ord] = sort(abs(coef1(active_peers)), 'descend');
                active_peers = active_peers(ord(1:min(max_peer_l3_per_equation, length(active_peers))));
            end
        end

        xobs2 = [];
        labels2 = {};

        for a = 1:length(active_inputs)
            ii = active_inputs(a);
            xobs2 = [xobs2, xobs(:,ii)];
            labels2{end+1} = mode_labels{ii};
        end

        for a = 1:length(active_peers)
            ii = active_peers(a);
            xobs2 = [xobs2, xobs(:,ii)];
            labels2{end+1} = mode_labels{ii};
        end

        xobs2 = [xobs2, xobs(:,target_col)];
        labels2{end+1} = mode_labels{target_col};
        target_col2 = length(labels2);

        base_for_quadratics = [active_inputs(:); active_peers(:); target_col];
        base_for_quadratics = unique(base_for_quadratics, 'stable');

        if length(base_for_quadratics) > max_quadratic_base_terms
            base_for_quadratics = base_for_quadratics(1:max_quadratic_base_terms);
        end

        if include_quadratic_terms
            for a = 1:length(base_for_quadratics)
                ii = base_for_quadratics(a);
                xobs2 = [xobs2, xobs(:,ii).^2];
                labels2{end+1} = sprintf('(%s)^2', mode_labels{ii});
            end
        end

        if include_cross_terms
            for a = 1:length(base_for_quadratics)
                for b = a+1:length(base_for_quadratics)
                    ii = base_for_quadratics(a);
                    jj = base_for_quadratics(b);

                    xobs2 = [xobs2, xobs(:,ii).*xobs(:,jj)];
                    labels2{end+1} = sprintf('(%s)*(%s)', mode_labels{ii}, mode_labels{jj});
                end
            end
        end

        [w2,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
            wsindy_ode_fun_capped(xobs2,tobs,weights, ...
            polys,trigs,custom_tags,custom_fcns, ...
            phi_class,max_d,tau,tauhat,K_frac,overlap_frac,relax_AG, ...
            scale_Theta,useGLS,lambda2,gamma,alpha_loss, ...
            overlap_frac_ag,pt_ag_fac,mt_ag_fac,run_sindy,useFD,smoothing_window);

        coef2 = w2(:, target_col2);

        coef_cap = 25;
        if any(abs(coef2) > coef_cap)
            fprintf('Clipping %s coefficients above %.1f in reduced model.\n', ...
                mode_labels{target_col}, coef_cap);
            coef2 = max(min(coef2, coef_cap), -coef_cap);
        end

        self_idx = find(strcmp(labels2, mode_labels{target_col}), 1);
        if ~isempty(self_idx) && coef2(self_idx) > 0
            fprintf('Removing positive self-feedback for %s: %.4e -> 0\n', ...
                mode_labels{target_col}, coef2(self_idx));
            coef2(self_idx) = 0;
        end

        coef2 = prune_equation_terms(coef2, labels2, mode_labels{target_col}, ...
            relative_term_tol, max_terms_per_equation, max_nonlinear_terms_per_equation);

        w_second{kk} = coef2;
        labels_second{kk} = labels2;
    end
end

function dydt = mrdmd_wsindy_inside_rhs(tt, y, w_second, labels_second, ...
    mode_labels, target_cols, xobs, tobs)

    dydt = zeros(length(target_cols),1);

    for kk = 1:length(target_cols)

        coef = w_second{kk};
        labels = labels_second{kk};

        theta = zeros(length(labels),1);

        for r = 1:length(labels)
            theta(r) = evaluate_inside_label(labels{r}, tt, y, ...
                mode_labels, target_cols, xobs, tobs);
        end

        dydt(kk) = coef(:).' * theta(:);
    end
end


function val = evaluate_inside_label(label, tt, y, mode_labels, target_cols, xobs, tobs)

    % self quadratic
    if startsWith(label, '(') && endsWith(label, ')^2')
        inner = extractBetween(label, "(", ")^2");
        base = evaluate_inside_label(inner{1}, tt, y, mode_labels, target_cols, xobs, tobs);
        val = base^2;
        return;
    end

    % cross term: (mode A)*(mode B)
    if contains(label, ')*(')
        parts = regexp(label, '^\((.*)\)\*\((.*)\)$', 'tokens');
        a = parts{1}{1};
        b = parts{1}{2};

        va = evaluate_inside_label(a, tt, y, mode_labels, target_cols, xobs, tobs);
        vb = evaluate_inside_label(b, tt, y, mode_labels, target_cols, xobs, tobs);

        val = va * vb;
        return;
    end

    % if label is a target L3 variable, use integrated y
    idx_target = find(strcmp(mode_labels(target_cols), label), 1);

    if ~isempty(idx_target)
        val = y(idx_target);
        return;
    end

    % otherwise, this is an input/ancestor mode.
    % Inside-bin recreation uses the observed xobs value at time tt.
    idx = find(strcmp(mode_labels, label), 1);

    if isempty(idx)
        error('Could not find label: %s', label);
    end

    val = interp1(tobs, xobs(:,idx), tt, 'linear', 'extrap');
end


function [lev, bin, mode_idx] = parse_mode_label(label)
nums = sscanf(label, 'L%d B%d Mode%d');

lev = nums(1);
bin = nums(2);
mode_idx = nums(3);
end

function [mean_l3_err, max_l3_err, mean_full_err, max_full_err] = evaluate_wsindy_parameter_set( ...
    list_w, list_b, list_modes, list_t_start, list_bin_widths, ...
    X, dt, m, n, input_levels, target_level, ...
    fit_start_idx, fit_end_idx, test_start_idx, test_end_idx, ...
    top_input_modes_per_level, top_target_modes, ...
    lambda1, lambda2, gamma, max_terms_per_equation, max_quadratic_base_terms)

    [w_second, labels_second, mode_labels, target_cols, mu_xobs, sigma_xobs] = run_mrdmd_wsindy( ...
        list_w, list_b, list_modes, list_t_start, list_bin_widths, ...
        dt, m, input_levels, target_level, fit_start_idx, fit_end_idx, ...
        top_input_modes_per_level, top_target_modes, lambda1, lambda2, gamma, ...
        max_terms_per_equation, max_quadratic_base_terms);
    
    test_steps = test_end_idx - test_start_idx + 1;
    recreate_steps = test_steps;
    
    [xobs_test, ~] = build_mrdmd_xobs_for_labels( ...
        list_w, list_b, list_t_start, list_bin_widths, ...
        dt, m, mode_labels, test_start_idx, test_end_idx);
    
    xobs_test = (xobs_test - mu_xobs) ./ sigma_xobs;
    tobs_test = ((test_start_idx:test_end_idx)' - 1) * dt;
    
    y_recreate = xobs_test(:, target_cols);
    
    for q = 1:(test_steps-1)
        y_now = xobs_test(q, target_cols).';
    
        dydt_now = mrdmd_wsindy_inside_rhs( ...
            tobs_test(q), y_now, ...
            w_second, labels_second, mode_labels, target_cols, ...
            xobs_test, tobs_test);
    
        y_recreate(q+1, :) = xobs_test(q, target_cols) + dt * dydt_now.';
    end
    
    y_recreate_phys = y_recreate .* sigma_xobs(target_cols) + mu_xobs(target_cols);
    
    X_l3_wsindy = zeros(n, recreate_steps, 'like', X);
    X_l3_true = zeros(n, recreate_steps, 'like', X);
    
    for k = 1:recreate_steps
        absolute_frame = test_start_idx + k - 1;
    
        for jj = 1:length(target_cols)
            mode_label = mode_labels{target_cols(jj)};
            [lev, bin, mode_idx] = parse_mode_label(mode_label);
    
            modes = list_modes{lev, bin};
            eigs_slow = list_w{lev, bin};
            b = list_b{lev, bin};
    
            t_start = list_t_start(lev, bin);
            rel_time = absolute_frame - t_start;
    
            X_l3_wsindy(:,k) = X_l3_wsindy(:,k) + ...
                modes(:, mode_idx) * y_recreate_phys(k,jj);
    
            X_l3_true(:,k) = X_l3_true(:,k) + ...
                modes(:, mode_idx) * (b(mode_idx) * eigs_slow(mode_idx)^rel_time);
        end
    end
    
    X_l3_true = real(X_l3_true);
    X_l3_wsindy = real(X_l3_wsindy);
    
    l3_recreate_error = zeros(1, recreate_steps);
    
    for k = 1:recreate_steps
        l3_recreate_error(k) = norm(X_l3_true(:,k) - X_l3_wsindy(:,k)) / ...
                                (norm(X_l3_true(:,k)) + eps);
    end
    
    X_inside_pred = zeros(n, recreate_steps, 'like', X);
    
    for lev = 1:2
        J = 2^(lev-1);
    
        for bin = 1:J
            if isempty(list_modes{lev,bin})
                continue;
            end
    
            t_start = list_t_start(lev,bin);
            bin_width = list_bin_widths(lev,bin);
            t_end = t_start + bin_width - 1;
    
            if t_start > test_end_idx || t_end < test_start_idx
                continue;
            end
    
            modes = list_modes{lev,bin};
            eigs_slow = list_w{lev,bin};
            b = list_b{lev,bin};
    
            for q = 1:recreate_steps
                absolute_frame = test_start_idx + q - 1;
    
                if absolute_frame < t_start || absolute_frame > t_end
                    continue;
                end
    
                rel_time = absolute_frame - t_start;
                time_powers = eigs_slow .^ rel_time;
    
                X_inside_pred(:,q) = X_inside_pred(:,q) + modes * (b .* time_powers);
            end
        end
    end
    
    X_inside_pred = real(X_inside_pred + X_l3_wsindy);
    X_inside_true = X(:, test_start_idx:test_end_idx);
    
    inside_full_error = zeros(1, recreate_steps);
    
    for k = 1:recreate_steps
        inside_full_error(k) = norm(X_inside_true(:,k) - X_inside_pred(:,k)) / ...
                                (norm(X_inside_true(:,k)) + eps);
    end
    
    mean_l3_err = mean(l3_recreate_error) * 100;
    max_l3_err = max(l3_recreate_error) * 100;
    mean_full_err = mean(inside_full_error) * 100;
    max_full_err = max(inside_full_error) * 100;

end

function coef = prune_equation_terms(coef, labels, target_label, relative_tol, max_terms, max_nonlinear_terms)
if isempty(coef) || all(coef == 0)
    return;
end

abs_coef = abs(coef);
max_coef = max(abs_coef);

if max_coef == 0
    return;
end

small_terms = abs_coef < relative_tol * max_coef;
coef(small_terms) = 0;

nonlinear = false(size(coef));
for ii = 1:length(labels)
    nonlinear(ii) = contains(labels{ii}, '^2') || contains(labels{ii}, ')*(');
end

active_nonlinear = find(coef ~= 0 & nonlinear);
if length(active_nonlinear) > max_nonlinear_terms
    [~, ord] = sort(abs(coef(active_nonlinear)), 'descend');
    drop_idx = active_nonlinear(ord(max_nonlinear_terms+1:end));
    coef(drop_idx) = 0;
end

active = find(coef ~= 0);
if length(active) > max_terms
    [~, ord] = sort(abs(coef(active)), 'descend');
    keep = active(ord(1:max_terms));

    target_self = find(strcmp(labels, target_label), 1);
    if ~isempty(target_self) && coef(target_self) ~= 0 && ~ismember(target_self, keep)
        weakest_keep = keep(end);
        keep(end) = target_self;
        fprintf('Keeping target self term for %s; dropping weaker term %.4e.\n', ...
            target_label, coef(weakest_keep));
    end

    drop = setdiff(active, keep);
    coef(drop) = 0;
end
end

function [xobs, tobs] = build_mrdmd_xobs_for_labels( ...
    list_w, list_b, list_t_start, list_bin_widths, ...
    dt, m, mode_labels, start_idx, end_idx)

m_interval = end_idx - start_idx + 1;
xobs = zeros(m_interval, length(mode_labels));

for c = 1:length(mode_labels)

    [lev, bin, mode_idx] = parse_mode_label(mode_labels{c});

    eigs_slow = list_w{lev, bin};
    b = list_b{lev, bin};
    bin_width = list_bin_widths(lev, bin);
    t_start = list_t_start(lev, bin);
    t_end = t_start + bin_width - 1;

    if start_idx > t_end || end_idx < t_start
        continue;
    end

    local_start = max(start_idx, t_start);
    local_end   = min(end_idx, t_end);

    rel_start = local_start - t_start;
    rel_end   = local_end - t_start;

    insert_start = local_start - start_idx + 1;
    insert_end   = local_end - start_idx + 1;

    a_local = b(mode_idx) * eigs_slow(mode_idx).^(rel_start:rel_end);
    xobs(insert_start:insert_end, c) = real(a_local(:));
end

tobs = (0:m_interval-1)' * dt;
end
