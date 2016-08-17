function [peaks troughs] = vt_phase_peak(cfg,indata)

% Computes the peak power of high frequency band depending on the peak and
% trough of a low frequency band
% Output are two structures, one with the trial x channel power of peaks
% and one with trial x channel power of troughs
%
% Input:
% cfg.lf = [Hz]; Low Frequency
% cfg.hf = [begin end]; High Frequency Range
% cfg.filttype = 'fir' or 'firls';
% cfg.toi = [begin end]; Time of Interest, default is whole interval
% cfg.per_win = Percent Window size around peak to average, default is +/-
% 5%
%
% (c) Julian Keil, 2014
% Version 1.0.: Recoded Peak/Trough-Script by Voytek with FieldTrip functions
%
% Based on PAC.m by
% Bradley Voytek
% Copyright (c) 2010
% University of California, Berkeley
% Helen Wills Neuroscience Institute

%% Set cfgs

low_frequency = cfg.lf; % define theta and alpha bands
high_frequency = cfg.hf; % define gamma band
filttype = cfg.filttype;
sampling_rate = indata.fsample;

%% Cut to desired latency
if isfield(cfg,'toi');
    tmpcfg = [];
    tmpcfg.toilim = cfg.toi;
    
    indata = ft_redefinetrial(tmpcfg,indata);
end

%% Define window around peaks and troughs to average gamma power
if isfield(cfg,'per_win')
    percent_window_size = cfg.per_win;
else
    percent_window_size = 0.05; % +/- 5 percent window around peaks and troughs
end

window_size = round(sampling_rate / low_frequency(1) .* percent_window_size);
clear percent_window_size

%% Filter the Data and get Phase
% Basics:
cfg=[];
cfg.hpfilter = 'yes';
cfg.hpfilttype = filttype;
cfg.lpfilter = 'yes';
cfg.lpfilttype = filttype;

cfg.hilbert = 'angle'; % OBACHT! Get the ANGLE!

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
cfg.hilbert = 'abs'; % OBACHT! Get the AMPLITUDE!
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

highdat = ft_preprocessing(cfg,indata); % Now we have the HIGHFREQ Signal

% Find peaks and troughs
for t = 1:length(lowdat.trial)
    fprintf('Searching High and Low on Trial %d of %d \n',t,length(lowdat.trial))
    for c = 1:size(lowdat.trial{t},1)
        % Find Troughs as local minima and peaks as zero crossings
        [lowMax, lowTrough, lowPeak] = findextrema(lowdat.trial{t}(c,:));
        % Throw out all *fake* zero crossings on the way from pi to -pi by
        % eliminating all lowPeak values that come right before a trough
        killvec = lowTrough -1;  
        lowPeak = setdiff(lowPeak,killvec);
        
        % Delete the first and last Trough if the windows would be outside
        % the EEG window
        [a] = find(lowTrough <= (window_size+1));
        if isempty(a)
            a=0;
        end
        lowTrough = lowTrough(max(a)+1:end); % Throw out the lows
        
        [a] = find(lowTrough >= length(lowdat.trial{t})-(window_size+1));
        if isempty(a)
            a=0;
        end
        lowTrough = lowTrough(1:end-max(a)); % Throw out the highs
        
        [a] = find(lowPeak <= (window_size+1));
        if isempty(a)
            a=0;
        end
        lowPeak = lowPeak(max(a)+1:end); % Throw out the lows
        
        [a] = find(lowPeak >= length(lowdat.trial{t})-(window_size+1));
        if isempty(a)
            a=0;
        end
        lowPeak = lowPeak(1:end-max(a)); % Throw out the highs
        
             
        % Now Average Gamma Power for Peaks and Troughs
        % Define Windows to average
        twin1 = lowTrough - window_size;
        twin2 = lowTrough + window_size;
        
        % Loop Windows
        for l = 1:length(twin1)
            tmp(l) = mean(highdat.trial{t}(c,[twin1(l):twin2(l)]));
        end
        
        % And Average
        troughPower(c) = mean(tmp);
            
        % Define Windows
         pwin1 = lowPeak - window_size;
         pwin2 = lowPeak + window_size;
        
         % Loop Windows
         for l = 1:length(pwin1)
            tmp(l) = mean(highdat.trial{t}(c,[pwin1(l):pwin2(l)]));
         end
         
         % Average
         peakPower(c) = mean(tmp);
            
         clear lowMax lowTrough lowPeak killvec twin* pwin*
    end
    pp{t} = peakPower;
    tp{t} = troughPower;
end

%% Put everything back into FT-Strcuture
peaks = indata;
peaks = rmfield(indata,'trial');
peaks.dimord = 'chan_time';
for t = 1:length(peaks.time)
    peaks.trial{t}(:,1) = pp{t};
    peaks.trial{t}(:,2) = pp{t};
    peaks.time{t} = [indata.time{t}(1) indata.time{t}(2)];
end

troughs = indata;
troughs = rmfield(indata,'trial');
troughs.dimord = 'chan_time';
for t = 1:length(troughs.time)
    troughs.trial{t}(:,1) = tp{t};
    troughs.trial{t}(:,2) = tp{t};
    troughs.time{t} = [indata.time{t}(1) indata.time{t}(2)];
end

end
