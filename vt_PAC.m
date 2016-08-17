function out = vt_PAC(cfg,indata)

% Computes Phase-Amplitude Coupling by estimating the phase of a low
% frequency oscillation and the phase of the envelope of a high frequency
% oscillation and subsequently computing PLV between both phases
%
% Input:
% cfg.lf = [begin end]; Low Frequency Range
% cfg.hf = [begin end]; High Frequency Range
% cfg.filttype = 'fir' or 'firls';
% cfg.toi = [begin end]; Time of Interest, default is whole interval
%
% (c) Julian Keil 2014
% Version 1.0.: Recoded PAC-Script by Voytek with FieldTrip functions
%
% Based on a function by:
% Bradley Voytek
% Copyright (c) 2010
% University of California, Berkeley
% Helen Wills Neuroscience Institute

%% Set cfgs

low_frequency = cfg.lf; % define theta and alpha bands
high_frequency = cfg.hf; % define gamma band
filttype = cfg.filttype;

%% Cut to desired latency
if isfield(cfg,'toi');
    tmpcfg = [];
    tmpcfg.toilim = cfg.toi;
    
    indata = ft_redefinetrial(tmpcfg,indata);
end

%% Filter the Data and get Phase
% Basics:
cfg=[];
cfg.hpfilter = 'yes';
cfg.hpfilttype = filttype;
cfg.lpfilter = 'yes';
cfg.lpfilttype = filttype;
cfg.hilbert = 'angle';

% LOW FREQS
% Formula (1) in Voytek 2010
disp('Low Frequency Filter');
cfg.hpfreq = low_frequency(1);
cfg.hpfiltord = 3*fix(indata.fsample/cfg.hpfreq);
    % Check if Filter Order is Odd
    if mod(cfg.hpfiltord,2) == 1
        disp('Filter Order is Odd, adding 1');
        cfg.hpfiltord = cfg.hpfiltord +1;
    end

cfg.lpfreq = low_frequency(2);
cfg.lpfiltord = 3*fix(indata.fsample/cfg.lpfreq);
    % Check if Filter Order is Odd
    if mod(cfg.lpfiltord,2) == 1
        disp('Filter Order is Odd, adding 1');
        cfg.lpfiltord = cfg.lpfiltord +1;
    end

lowdat = ft_preprocessing(cfg,indata);

% HIGH FREQS -> Filter and then Plug into Filter again (Filterception)
% Formula (1) in Voytek 2010
disp('High Frequency Filter');
cfg.hpfreq = high_frequency(1);
cfg.hpfiltord = 3*fix(indata.fsample/cfg.hpfreq);
    if mod(cfg.hpfiltord,2) == 1
        % Check if Filter Order is Odd
        disp('Filter Order is Odd, adding 1');
        cfg.hpfiltord = cfg.hpfiltord +1;
    end
    
cfg.lpfreq = high_frequency(2);
cfg.lpfiltord = 3*fix(indata.fsample/cfg.lpfreq);
    if mod(cfg.lpfiltord,2) == 1
        % Check if Filter Order is Odd
        disp('Filter Order is Odd, adding 1');
        cfg.lpfiltord = cfg.lpfiltord +1;
    end

tmp = ft_preprocessing(cfg,indata); % Now we have the HIGHFREQ Signal

% Plug HIGHFREQ back into LOWFREQ Filter
% Formula (2) in Voytek 2010
disp('High Frequency Envelope');
cfg.hpfreq = low_frequency(1);
cfg.hpfiltord = 3*fix(indata.fsample/cfg.hpfreq);
    % Check if Filter Order is Odd
    if mod(cfg.hpfiltord,2) == 1
        disp('Filter Order is Odd, adding 1');
        cfg.hpfiltord = cfg.hpfiltord +1;
    end

cfg.lpfreq = low_frequency(2);
cfg.lpfiltord = 3*fix(indata.fsample/cfg.lpfreq);
    % Check if Filter Order is Odd
    if mod(cfg.lpfiltord,2) == 1
        disp('Filter Order is Odd, adding 1');
        cfg.lpfiltord = cfg.lpfiltord +1;
    end

highdat = ft_preprocessing(cfg,tmp);

%% Calculate PAC

% Loop Trials
for t=1:length(indata.trial);
    %fprintf('Working on trial %d of %d \n',t,length(indata.trial))
    
    % Initialize PAC variable
    pac{t} = zeros([size(indata.trial{t}, 1)]);
    
    % Formula (3) in Voytek 2010
    pac{t} = abs(sum(exp(1i * (lowdat.trial{t} - highdat.trial{t})),2))/size(highdat.trial{t},2);
    
end

%% Write back into FT-Strcuture

out = indata;

for t = 1:length(out.time)
    out.trial{t} = [pac{t} pac{t}];
    out.time{t} = [indata.time{t}(1) indata.time{t}(2)];
end