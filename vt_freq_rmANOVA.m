function stats = vt_freq_rmANOVA(cfg,varargin)
%% ANOVA TEST
% Compute 2-Factor Within Subject Condition and Interactions for all EEG-Sensors
%
% cfg.channel = channel indices/names
% cfg.avgoverchan = Average over Channels, or compute channel-wise ANOVA
% cfg.avgoverfreq = Average over Frequency
% cfg.freq;
% cfg.IV1 = number of Levels in Factor 1
% cfg.IV2 = number of Levels in Factor 2
% cfg.parameter
% cfg.neighbours = neighbourhood structure for spatial correction
% cfg.correctm = 'neighbours' or 'fdr';
% cfg.minnb = minimum number of neighbours, default = 1;
%
% USE AS:
%cfg=[];
% cfg.frequency = [8 12];
% cfg.avgoverfreq = 'yes';
% cfg.nIV1 = 2;
% cfg.nIV2 = 2;
% cfg.parameter = 'powspctrm';
% cfg.alpha = .001;
% 
% stats = vt_freq_rmANOVA(cfg,fft_vf{:},fft_tf{:},fft_vv{:},fft_tv{:});
%
% (c) Julian Keil 16.04.2014
% Version 1: Implemented Function
% Version 2: Changed name to Freq
% Version 3: Corrected Error in Averaging
% Version 4: TimeFrequency-Version
% Version 5: Updated channel-Selection Options (04.11.2014)
% Version 5.1.: Updated Channel Selection for Neighbour-Correction (12.12.2014
% Version 5: Major Bugfix in Neighbour-Correction (04.11.2015) 
% Version 5.1.: Minor Bugfix in Channel options and corection (30.05.2016)
% Version 6.: Bugfix in the preparation for the Anova-Matrix. Now supports
% nXm-Anovas (15.06.2017)

%% Make Freq to Time-Freq in case of FFT
if strcmpi(varargin{1}.dimord,'chan_freq')
    for v = 1:length(varargin)
        %varargin{v}.powspctrm = [varargin{v}.powspctrm ];
        varargin{1}.time = 1;
    end
end
%% Set CFGs

nIV1 = cfg.nIV1; % Factors on Level1
nIV2 = cfg.nIV2; % Factors on Level2
nsub = numel(varargin)/(nIV1*nIV2); % How many subjects in total?

% Set Alpha Level
if isfield(cfg,'alpha')
    alpha = cfg.alpha;
else
    alpha = .05;
end

% Set Correction 
if isfield(cfg,'correctm') % Check for cfg
    correctm = cfg.correctm; % set Method
    if strcmpi(correctm,'fdr') 
        fprintf('Using FDR Correction \n')
    elseif strcmpi(correctm,'neighbours') && isfield(cfg,'neighbours')
        neigh = cfg.neighbours; % Use Predefined Neighbours
        fprintf('Using Neighbour Correction \n')
    else
        fprintf('Please use either FDR or provide Neighbourhood Structure \n');
        return
    end
    % Set Minimum Neighbours
    if strcmpi(correctm,'neighbours') && isfield(cfg,'minnb')
        minnb = cfg.minnb;
    else
        minnb = 1;
    end
    
else
    correctm = 'none';
end

% Set Time Options
if isfield(cfg,'latency');
    if isfield(cfg,'avgovertime');
        if strcmpi(cfg.avgovertime,'yes');
            timerange = [nearest(varargin{1}.time,cfg.latency(1)) nearest(varargin{1}.time,cfg.latency(end))];
        end
    else
        timerange = [nearest(varargin{1}.time,cfg.latency(1)) nearest(varargin{1}.time,cfg.latency(end))]; 
        latency = varargin{1}.time(timerange(1):timerange(end));
    end
else       
    ntime = length(varargin{1}.time); % How many Time Steps
    timerange = 1:ntime;   
    latency = varargin{1}.time(timerange(1):timerange(end));
end

% Set Frequency Options
if isfield(cfg,'frequency');
    if isfield(cfg,'avgoverfreq');
        if strcmpi(cfg.avgoverfreq,'yes');
            freqrange = [nearest(varargin{1}.freq,cfg.frequency(1)) nearest(varargin{1}.freq,cfg.frequency(end))];
        end
    else
    freqrange = [nearest(varargin{1}.freq,cfg.frequency(1)) nearest(varargin{1}.freq,cfg.frequency(end))];
    frequency = varargin{1}.freq(freqrange(1):freqrange(end));
    end
else       
    nfreq = length(varargin{1}.freq); % How many Freq Steps
    freqrange = 1:nfreq;    
    frequency = varargin{1}.freq(freqrange(1):freqrange(end));
end

if isfield(cfg,'parameter')
    param = cfg.parameter; % What do we want to work on
else
    param = 'powspctrm';
end


% Set Channel Options
if isfield(cfg,'channel');
    if strcmpi(cfg.channel,'all')
        nchan = size(varargin{1}.(param),1);
        inchan = 1:nchan;
    else
        nchan = length(cfg.channel);
        if iscell(cfg.channel)
            for c = 1:length(cfg.channel)
            inchan(c) = find(strcmpi(cfg.channel{c},varargin{1}.label)==1);
            end
        else
%             tmp =  cellfun(@num2str,num2cell(cfg.channel),'UniformOutput',false);
%             for c = 1:length(tmp)
%                 inchan(c) = find(strcmpi(tmp{c},varargin{1}.label)==1);
%             end
            inchan=cfg.channel;
        end
        
    end
else
    nchan = size(varargin{1}.(param),1);
    inchan = 1:nchan;
end
    
%% Prepare Data MAtrix
% Empty Matrix for Data
Mat = zeros(nsub*nIV1*nIV2,4); % Subjects X Conditons

% IV1
IV1 = [];
for g=1:nIV1
    tmp = repmat((ones(nsub,1)*g),nIV2,1);
    IV1 = [IV1; tmp ];
end
Mat(:,2) = IV1; % IV1

% IV2
IV2 = [];
for c=1:nIV2
    tmp = ones(nsub,1)*c;
    IV2 = [IV2; tmp ];
end
Mat(:,3) = repmat(IV2,nIV1,1); % IV2

% Subjects
Subj = [];
for s=1:nIV1*nIV2;
    tmp = [1:nsub]';
    Subj = [Subj; tmp];
end
Mat(:,4) =  Subj; % Subject X Levels
%% Plug in Data
% Look for the right element and put it into the first column of Mat
tic
% Should we average over channels?
if isfield(cfg,'avgoverchan')
    if strcmpi(cfg.avgoverchan,'yes')
        for v = 1:length(varargin)
            varargin{v}.(param) = mean(varargin{v}.(param)(inchan,:,:));
        end
        nchan = 1; % Update Number of Channels to 1!
        inchan = 1; % Update Channel Input to 1
    end
end

% Should we average over time
if isfield(cfg,'avgovertime')
    if strcmpi(cfg.avgovertime,'yes')
        for v = 1:length(varargin)
            varargin{v}.(param) = mean(varargin{v}.(param)(:,:,timerange(1):timerange(end)),3);
            varargin{1}.time = 1;
        end
        latency = 1;
    end
end

% Should we average over frequency
if isfield(cfg,'avgoverfreq')
    if strcmpi(cfg.avgoverfreq,'yes')
        for v = 1:length(varargin)
            varargin{v}.(param) = mean(varargin{v}.(param)(:,freqrange(1):freqrange(end),:),2);
            varargin{1}.freq = 1;
        end
        frequency = 1;
    end
end
%% Let's go!

% Pre-allocate Variables
tmp_g = NaN(size(Mat(:,1)));
means = NaN([nchan,length(frequency),length(latency),nIV1,nIV2]);
sems = means;
            
for n = 1:nchan % For Channels
    fprintf('Statistical Stuff at Channel %i of %i \n', n, nchan)
    for t = 1:length(latency);
        fprintf('Statistical Stuff at Timepoint %i of %i \n', t, length(latency))
        for f = 1:length(frequency);
            
            fi = nearest(varargin{1}.freq,frequency(f));
            ti = nearest(varargin{1}.time,latency(t));
            % Select the relevant elements
            for c=1:nIV1*nIV2 % Loop Conditions
                for s=1:nsub % loop Subjects varagin
                    sid = sub2ind([nsub, nIV1*nIV2],s,c); % Subject Identifier, Subjects are order 1...n;1...n.. in varargin 
                    tmp_g(sid,1) = varargin{sid}.(param)(inchan(n),fi,ti);  % take the relevant element from the input
                end % for Conditions  
            end % for Subjects

            Mat(:,1) = tmp_g; % MAtrix for ANOVA

            %% Descriptive Stats       
            for g = 1:nIV1
                for c = 1:nIV2 % 
                    % Output is ChannelXFreqXIV1XIV2
                    tmp = mean(Mat((Mat(:,3)==c) & (Mat(:,2)==g),:));
                    means(n,f,t,g,c) = tmp(1); 
                    tmp_z = std(Mat((Mat(:,3)==c) & (Mat(:,2)==g) ,:));
                    tmp_n = sqrt(length(Mat(Mat(:,3)==c & Mat(:,2)==g)));
                    sems(n,f,t,g,c) = tmp_z(1)/tmp_n(1);
                end
            end
            
            %% Compute the Anova
            % Blantantly used the existing function but commented the output out
            %  Trujillo-Ort
            out = vt_RMAOV2(Mat,alpha); % Hacked Anova Function
            
             %% Put into little boxes

            stats.statIV1(n,f,t) = out.F1;
            stats.probIV1(n,f,t) = out.P1;
            stats.dfIV1(:,n,f) = [out.v1, out.v2];
            stats.statIV2(n,f,t) = out.F2;
            stats.probIV2(n,f,t) = out.P2;
            stats.dfIV2(:,n,f) = [out.v3, out.v2];
            stats.statint(n,f,t) = out.F3;
            stats.probint(n,f,t) = out.P3;
            stats.dfint(:,n,f) = [out.v5, out.v6];

        end
    end
end

stats.means = means;
stats.sems = sems;
stats.dimord = 'chan_freq_time';

%% Correction
switch correctm
    case{'none'}
        fprintf('Not performing correction for multiple comparisons \n Mask is uncorrected! \n')
            stats.maskIV1 = stats.probIV1<alpha;
            stats.maskIV2 = stats.probIV2<alpha;
            stats.maskint = stats.probint<alpha;
    case{'fdr'} % Using the FDR method 
        %  Arnaud Delorme
        fprintf('Using FDR correction for multiple comparisons \n')
            [stats.probIV1_fdr, stats.maskIV1] = vt_fdr(stats.probIV1, alpha);
            [stats.probIV2_fdr, stats.maskIV2] = vt_fdr(stats.probIV2, alpha);
            [stats.probint_fdr, stats.maskint] = vt_fdr(stats.probint, alpha);
    case{'neighbours'}
        fprintf('Using Neighbours to correct for for multiple comparisons \n')
            stats.maskIV1 = zeros(size(stats.statIV1));
            stats.maskIV2 = zeros(size(stats.statIV2));
            stats.maskint = zeros(size(stats.statint));  
            
            % Get Input Labels
            if iscell(varargin{1}.label)
                lab = varargin{1}.label;
            else
                lab = cellfun(@str2num,varargin{1}.label,'UniformOutput',false); 
            end
       
        for c = 1:length(lab)
            for f = 1:length(frequency)
                for t = 1:length(latency)
                    
                if stats.probIV1(c,f,t) < alpha % is smaller than alpha?
                    if iscell(neigh(c).neighblabel)
                        nb = neigh(c).neighblabel;
                        % Get Neighbour Indices
                        for cn = 1:length(nb)
                            a(cn) = find(strcmpi(nb(cn),lab));
                        end
                    else
                        nb =  cellfun(@str2num,neigh(c).neighblabel,'UniformOutput',false);
                        % Get Neighbour Indices
                        for cn = 1:length(nb)
                            a(cn) = find(nb(cn) == lab);
                        end
                    end
                    % Sort the Indexed Channels
                    nb_alpha = sort(stats.probIV1(a,f,t)); 
                    if nb_alpha(minnb) < alpha % is the nth neighbour also smaller alpha?
                        stats.maskIV1(c,f,t) = 1;
                    else
                        stats.maskIV1(c,f,t) = 0;
                    end % If Neighbour-Alpha
                end % If Group

                clear a nb nb_alpha

                if stats.probIV2(c,f,t) < alpha % Condition
                    if iscell(neigh(c).neighblabel)
                        nb = neigh(c).neighblabel;
                        % Get Neighbour Indices
                        for cn = 1:length(nb)
                            a(cn) = find(strcmpi(nb(cn),lab));
                        end
                    else
                        nb =  cellfun(@str2num,neigh(c).neighblabel,'UniformOutput',false);
                        % Get Neighbour Indices
                        for cn = 1:length(nb)
                            a(cn) = find(nb(cn) == lab);
                        end
                    end
                    % Sort the Indexed Channels
                    nb_alpha = sort(stats.probIV2(a,f,t)); 
                    if nb_alpha(minnb) < alpha
                        stats.maskIV2(c,f,t) = 1;
                    else
                        stats.maskIV2(c,f,t) = 0;
                    end % If Neighbour-Alpha
                end % If Cond

                clear a nb nb_alpha

                if stats.probint(c,f,t) < alpha % Interaction
                    if iscell(neigh(c).neighblabel)
                        nb = neigh(c).neighblabel;
                        % Get Neighbour Indices
                        for cn = 1:length(nb)
                            a(cn) = find(strcmpi(nb(cn),lab));
                        end
                    else
                        nb =  cellfun(@str2num,neigh(c).neighblabel,'UniformOutput',false);
                        % Get Neighbour Indices
                        for cn = 1:length(nb)
                            a(cn) = find(nb(cn) == lab);
                        end
                    end
                    % Sort the Indexed Channels
                    nb_alpha = sort(stats.probint(a,f,t)); 
                    if nb_alpha(minnb) < alpha
                        stats.maskint(c,f,t) = 1;
                    else
                        stats.maskint(c,f,t) = 0;
                    end % If Neighbour-Alpha
                end % If Interaction

                clear a nb nb_alpha

                end % For Time
            end % For Freq
        end % For Chan
end % Switch
toc

%% add FT-Stuff
stats.freq = frequency;
stats.time = latency;
stats.label = varargin{1}.label(inchan);

