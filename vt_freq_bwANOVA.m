function stats = vt_freq_bwANOVA(cfg,varargin)
%% ANOVA TEST
% Compute Group, Condition and Interactions for all EEG-Sensors
% Main Effects for Group are identical to ft_freqstatistics either for
% single frequencies or averge over frequencies
%
% cfg.channel = channel indices (NOT NAMES!)
% cfg.avgoverchan = Average over Channels, or compute channel-wise ANOVA
% cfg.avgovertime = Average over Interval or Compute Pointwise
% cfg.avgoverfreq = Average over Freuqencies oder Compute Freqwise
% cfg.latency = Time window to compute ANOVA;
% cfg.frequenc = Frequency Window to compute ANOVA
% cfg.ngr = number of groups
% cfg.sgr = size of groups as vector
% cfg.ncond = number of conditions
% cfg.parameter
% cfg.neighbours = neighbourhood structure for spatial correction
% cfg.correctm = 'neighbours' or 'fdr';
% cfg.minnb = minimum number of neighbours, default = 1;
%
% USE AS:
% cfg = [];
% cfg.channel = 'all';
% cfg.avgoverchan = 'no';
% cfg.latency = [.1 .15];
% cfg.frequency = []
% cfg.ngr = 2;
% cfg.sgr = [13 20];
% cfg.ncond = 4;
% cfg.parameter = 'avg';
% 
% stats = vt_freq_bwANOVA(cfg,p_ill{:},...
%     p_noill{:},...
%     c_ill{:},..
%     c_noill{:});
%
% (c) Julian Keil 13.12.2013
% Version 1.1.: 12.02.2014 - Added Channel Option
% Version 2.: 08.05.2014 - Modified to Timelock Data
% Version 3.: 12.05.2014 - Added Serial Anova (Timepoint-Wise) Option
% Version 4.: 23.05.2014 - Adapted time-function for time-freq
% Version 5.: 30.10.2014 - Major Bugfix in varargin sorting
% Version 6.: 19.05.2015 - Bug in Freqrange, Text output and Neighbour correction fixed
% Version 6.1.: 30.05.2016 - Minor Bugs in Channel selection and Correction
%% Set CFGs

latency=cfg.latency; % Latency
frequency=cfg.frequency; % Frequency

ngr = cfg.ngr; % Number of Groups
sgr = cfg.sgr; % Size of Groups
nsub = sum(sgr); % How many subjects in total?
ncond = cfg.ncond; % Number of Conditions


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
            ntime = 1;
            timerange = [nearest(varargin{1}.time,cfg.latency(1)) nearest(varargin{1}.time,cfg.latency(end))];
        end
    else
    ntime = length(cfg.latency); % How many Time Steps
    timerange = [nearest(varargin{1}.time,cfg.latency(1)) nearest(varargin{1}.time,cfg.latency(end))]; 
    end
else       
    ntime = length(varargin{1}.time); % How many Time Steps
    timerange = 1:ntime;    
end

% Set Frequency Options
if isfield(cfg,'frequency');
    if isfield(cfg,'avgoverfreq');
        if strcmpi(cfg.avgoverfreq,'yes');
            nfreq = 1;
            freqrange = [nearest(varargin{1}.freq,cfg.frequency(1)) nearest(varargin{1}.freq,cfg.frequency(end))];
        end
    else
    nfreq = length(cfg.frequency); % How many Time Steps
    freqrange = [nearest(varargin{1}.freq,cfg.frequency(1)) nearest(varargin{1}.freq,cfg.frequency(end))]; 
    end
else       
    nfreq = length(varargin{1}.freq); % How many Freq Steps
    freqrange = 1:nfreq;    
end

if isfield(cfg,'parameter')
    param = cfg.parameter; % What do we want to work on
else
    param = 'powspctrm';
end


ncond = cfg.ncond; % How many Freqs
param = cfg.parameter; % What do we want to work on

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
Mat = zeros(nsub*ncond,4); % Subjects X Conditons

% Groups
Gr = [];
for g=1:ngr
    tmp = ones(sgr(g)*ncond,1)*g;
    Gr = [Gr; tmp ];
end
Mat(:,2) = Gr; % Group

% Conditions -> That can be your Frequencies
Con = [];
for c=1:ncond
    tmp = ones(nsub,1)*c;
    Con = [Con; tmp];
end
tmp = 1:ncond;
Con = repmat(tmp',nsub,1);
Mat(:,3) = Con; % Condition

% Subjects
Subj = [];
for s=1:nsub;
    tmp = repmat(s,ncond,1);
    Subj = [Subj; tmp];
end
Mat(:,4) =  Subj; % 20 Subject X 3 Cond

%% Plug in Data
% Look for the right element and put it into the first column of Mat
tic
% Should we average over channels?
if isfield(cfg,'avgoverchan')
    if strcmpi(cfg.avgoverchan,'yes')
        for v = 1:length(varargin)
            varargin{v}.(param) = nanmean(varargin{v}.(param)(inchan,:,:));
        end
        nchan = 1; % Update Number of Channels to 1!
        inchan = 1; % Update Channel Input to 1
    end
end

% Should we average over time
if isfield(cfg,'avgovertime')
    if strcmpi(cfg.avgovertime,'yes')
        for v = 1:length(varargin)
            varargin{v}.(param) = nanmean(varargin{v}.(param)(:,:,timerange(1):timerange(end)),3);
            varargin{1}.time = 1;
        end
        latency = 1;
    end
end

% Should we average over frequency
if isfield(cfg,'avgoverfreq')
    if strcmpi(cfg.avgoverfreq,'yes')
        for v = 1:length(varargin)
            varargin{v}.(param) = nanmean(varargin{v}.(param)(:,freqrange(1):freqrange(end),:),2);
            varargin{1}.freq = 1;
        end
        frequency = 1;
    end
end
%% Let's go!
for n = 1:nchan % For Channels
    fprintf('Statistical Stuff at Channel %i of %i \n', n, nchan)
    for t = 1:length(latency);
        fprintf('Statistical Stuff at Timepoint %i of %i \n', t, length(latency))
        for f = 1:length(frequency);
        %fprintf('Statistical Stuff at Frequency %i of %i \n', f, length(frequency))

        tmpind = [];
        for g = 1:ngr
            tmp = repmat(sgr(g),1,sgr(g)*ncond);
            tmpind = [tmpind,tmp];
        end

        varind = Mat(:,4);
         % First Half
         for i = 2:2:length(varind)/2
            varind(i) = Mat(i,4)+tmpind(i);
         end
         
         % Second Half
         % First add length of group to all
         for i = ((length(varind)/2)+1):1:length(varind)
             varind(i) = Mat(i,4)+tmpind(i);
         end
         % Then add another length of group to every second item
         for i = ((length(varind)/2)+2):2:length(varind)
             varind(i) = varind(i)+tmpind(i);
         end
         
         
        ti = nearest(varargin{1}.time,latency(t));
        fi = nearest(varargin{1}.freq,frequency(f));
        for s = 1:length(varind)
            Mat(s,1) = varargin{varind(s)}.(param)(inchan(n),fi,ti);
        end
        dat = Mat(:,1);

        %% Descriptive Stats
        for g = 1:ngr
            for c = 1:ncond % Freq in this case!
                % Output is ChannelXGroupXFreq
                means(n,f,t,g,c) = squeeze(mean(dat((Mat(:,3)==c) & (Mat(:,2)==g) ,:)));
                sems(n,f,t,g,c) = squeeze(std(dat((Mat(:,3)==c) & (Mat(:,2)==g) ,:)))/(sqrt(length(dat((Mat(:,3)==c) & (Mat(:,2)==g) ,:))));
            end
        end

        %% Compute the Anova
        % Blantantly used the existing function but commented the output out
        %  Trujillo-Ortiz
        out = vt_BWAOV2(Mat,alpha); % Hacked Anova Function

         %% Put into little boxes

        stats.statgroup(n,f,t,:) = out.F1;
        stats.probgroup(n,f,t,:) = out.P1;
        stats.dfgroup(:,n) = [out.v1, out.v2];
        stats.statcond(n,f,t,:) = out.F2;
        stats.probcond(n,f,t,:) = out.P2;
        stats.dfcond(:,n) = [out.v3, out.v2];
        stats.statint(n,f,t,:) = out.F3;
        stats.probint(n,f,t,:) = out.P3;
        stats.dfint(:,n) = [out.v4, out.v5];
        end % Freq
    end % Time
end % Chan

stats.means = means;
stats.sems = sems;

%% Correction
switch correctm
    case{'none'}
        fprintf('Not performing correction for multiple comparisons \n Mask is uncorrected! \n')
            stats.maskgroup = stats.probgroup<alpha;
            stats.maskcond = stats.probcond<alpha;
            stats.maskint = stats.probint<alpha;
    case{'fdr'} % Using the FDR method 
        %  Arnaud Delorme
        fprintf('Using FDR correction for multiple comparisons \n')
            [stats.probgroup_fdr, stats.maskgroup] = vt_fdr(stats.probgroup, alpha);
            [stats.probcond_fdr, stats.maskcond] = vt_fdr(stats.probcond, alpha);
            [stats.probint_fdr, stats.maskint] = vt_fdr(stats.probint, alpha);
    case{'neighbours'}
        fprintf('Using Neighbours to correct for for multiple comparisons \n')
        stats.maskgroup = zeros(size(stats.statgroup));
        stats.maskcond = zeros(size(stats.statcond));
        stats.maskint = zeros(size(stats.statint));
        lab = cellfun(@str2num,varargin{1}.label); % Get Input Labels
        for c = 1:size(stats.statgroup,1)
            for f = 1:size(stats.statgroup,2)
                for t = 1:size(stats.statgroup,3)
               % fprintf('Checking Channel %d \n',c)
                if stats.probgroup(c,f,t,:) < alpha % Group
                    nb =  cellfun(@str2num,neigh(c).neighblabel);
                    % Get Neighbour Indices
                    for cn = 1:length(nb)
                        a(cn) = find(nb(cn) == lab);
                    end
                    % Sort the Indexed Channels
                    nb_alpha = sort(stats.probgroup(a,:)); 
                    if nb_alpha(minnb) < alpha
                        stats.maskgroup(c,f,t,:) = 1;
                    else
                        stats.maskgroup(c,f,t,:) = 0;
                    end % If Neighbour-Alpha
                end % If Group

                clear a nb nb_alpha

                if stats.probcond(c,f,t,:) < alpha % Condition
                    nb =  cellfun(@str2num,neigh(c).neighblabel);
                    % Get Neighbour Indices
                    for cn = 1:length(nb)
                        a(cn) = find(nb(cn) == lab);
                    end
                    % Sort the Indexed Channels
                    nb_alpha = sort(stats.probcond(a,:)); 
                    if nb_alpha(minnb) < alpha
                        stats.maskcond(c,f,t,:) = 1;
                    else
                        stats.maskcond(c,f,t,:) = 0;
                    end % If Neighbour-Alpha
                end % If Cond

                clear a nb nb_alpha

                if stats.probint(c,f,t,:) < alpha % Interaction
                    nb =  cellfun(@str2num,neigh(c).neighblabel);
                    % Get Neighbour Indices
                    for cn = 1:length(nb)
                        a(cn) = find(nb(cn) == lab);
                    end
                    % Sort the Indexed Channels
                    nb_alpha = sort(stats.probint(a,:)); 
                    if nb_alpha(minnb) < alpha
                        stats.maskint(c,f,t,:) = 1;
                    else
                        stats.maskint(c,f,t,:) = 0;
                    end % If Neighbour-Alpha
                end % If Interaction

                clear a nb nb_alpha
                end % For Time
            end % For Freq
        end % For Channel

end % Switch
toc

%% add FT-Stuff
stats.dimord = varargin{1}.dimord;
stats.time = latency;
stats.freq = frequency;
stats.label = varargin{1}.label(inchan);

