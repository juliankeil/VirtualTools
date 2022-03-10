function [data_in] = vt_autoreject(cfg,data_in)

% This function first removes bad channels (across trials) and then bad trials (across the remaining channels)
% The input requires raw data.
% 
% The following options are available
% cfg.keeptrials = 'yes' or 'no'; % Should bad trials be removed?
% cfg.keepchannels = 'yes' or 'no'; % Should bad channels be removed?
% cfg.threshold = numeric; % The outlier threshold (e.g. 2 means 2
% standardeviations above the mean for each metric)
% cfg.peak = numeric; % The absolute peak +/- Amplitude
%
% Julian Keil, 2022
% Ver 1.: 10.03.2022: First implementation

%% 0. Set the cfgs
if exist('cfg')
    % Trials
    if isfield(cfg,'keeptrials')
        if strcmpi(cfg.keeptrials,'no')
            trial_flag = 1;
        elseif strcmpi(cfg.keeptrials,'yes')
            meth_flag = 2;
        end
    else
        trial_flag = 1; % Default: Don't keep trials
    end
    %Channels
    if isfield(cfg,'keepchannels')
        if strcmpi(cfg.keepchannels,'no')
            chan_flag = 1;
        elseif strcmpi(cfg.keepchannels,'yes')
            chan_flag = 2;
        end
    else
        chan_flag = 1; % Default: Don't keep channels
    end
    % Thresholds
    if isfield(cfg,'threshold')
        thresh = cfg.threshold;
    else
        thresh = 2.5; % Default is 2.5 STD above the mean
    end
    % Peaks
    if isfield(cfg,'peak')
        peak = cfg.peak;
    else
        peak = 1000; % Default is 1000mV
    end
end % exist

clear goodchannels goodtrials
%% 1. First. remove the Channels 
% We'll restrict this to outliers based on 2.5 STD above the mean
% First, compute the std across trials for each channel
% Then, compute the mean and std of the stds across channels
if chan_flag == 1
    % 1.1. STD across trials
    tmp_c = zeros(1,length(data_in.trial)*length(data_in.trial{1})); % Build an empty vector to save time
    for c = 1:length(data_in.label) % Loop channels
        ti = 1; % Index of trials
        for t = 1:length(data_in.trial) % Loop trials
            tmp_c(ti:ti+length(data_in.trial{t})-1) = data_in.trial{t}(c,:); % Stitch all trials together
            ti = ti+length(data_in.trial{t}); % Move the trial index
        end
        std_c(c) = std(tmp_c); % And compute the std across trials
    end

    % 1.2. Compute Threshold
    m_std = mean(std_c); % Mean of STD across trials
    s_std = std(std_c); % STD of STD across trials

    % 1.3. Find good channels
    goodchannels = std_c(:) <= m_std + (thresh*s_std); % Mean + 2.5*STD            
    fprintf('\n Good channels: %d \n\n',length(goodchannels)); 

    % 3.1.1.4. Keep only good channels
    cfg = [];
    cfg.channel = data_in.label(logical(goodchannels));

    data_in = ft_selectdata(cfg,data_in);
end
        
%% Second, clean the trials
% For each trial, we'll check the kurtosis, z-value and mean across
% channels
% First, we'll compute the indices across channels for each trial
% Then, we'll remove trials based on thresholds
if trial_flag == 1
    % 3.1.2.1. Indices across channels
    tmp_t = zeros(1,length(data_in.label)*length(data_in.trial{1})); % Build an empty vector to save time
    for t = 1:length(data_in.trial) % Loop trials
        ci = 1; % Index of channels
        for c = 1:length(data_in.label) % Loop channels
            tmp_t(ci:ci+length(data_in.trial{t})-1) = data_in.trial{t}(c,:); % Stitch all trials together
            ci = ci+length(data_in.trial{t}); % Move the trial index
        end
        std_t(t) = std(tmp_t); % And compute the std across trials
        kurt_t(t) = kurtosis(tmp_t,[],2); % Kurtosis
        [zsc_t, mu_t, sig_t] = zscore(tmp_t,0,2); % Zscore
        zval_t(t) = max(abs(zsc_t)); % Max Zscore
        max_t(t) = max(tmp_t); % Peak Max
        min_t(t) = min(tmp_t); % Peak Min
    end

    % 3.1.2.2. Define thresholds: 2.5 STD above the mean
    z_thresh = mean(zval_t) + thresh*std(zval_t);
    k_thresh = mean(kurt_t) + thresh*std(kurt_t);
    v_thresh = mean(std_t) + thresh*std(std_t);

    % 3.1.2.3. Apply thresholds
    clear outs
    for t = 1:length(data_in.trial)
        outs(t) = any(zval_t(t) >= z_thresh ...
                    | kurt_t(t) >= k_thresh ...
                    | std_t(t) >= v_thresh ...
                    | max_t(t) == peak ...
                    | min_t(t) == -peak);
    end

    goodtrials = find(outs == 0); % Select what's left
    fprintf('\n Good trials: %d \n\n',length(goodtrials)); 
    % 3.1.2.4. Keep only good trials
    cfg = [];
    cfg.trials = goodtrials;

    data_in = ft_selectdata(cfg,data_in); % data_psc is the, preprocessed, selected and clean data
end