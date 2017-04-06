function stats = vt_time_bwANOVA(cfg,varargin)
%% ANOVA TEST
% Compute Group, Condition and Interactions for all EEG-Sensors
% Main Effects for Group are identical to ft_freqstatistics either for
% single frequencies or averge over frequencies
%
% cfg.channel = channel indices (NOT NAMES!)
% cfg.avgoverchan = Average over Channels, or compute channel-wise ANOVA
% cfg.latency;
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
% cfg.ngr = 2;
% cfg.sgr = [13 20];
% cfg.ncond = 4;
% cfg.parameter = 'avg';
% 
% stats = vt_time_bwANOVA(cfg,p_ill{:},...
%     p_noill{:},...
%     c_ill{:},..
%     c_noill{:});
%
% (c) Julian Keil 13.12.2013
% Version 1.1.: 12.02.2014 - Added Channel Option
% Version 2.: 08.05.2014 - Modified to Timelock Data
% Version 3.: 12.05.2014 - Added Serial Anova (Timepoint-Wise) Option
% Version 4.: Major BugFix in Varargin Sorting
% Version 4.1.: 30.05.2016 Minor Bugfixes in Correction and Channel Option
% Version 4.2.: 15.03.2017 - Minor Bugs in CFGs
% Version 5.: 15.03.2017 Another Major BugFix in Varargin Sorting 
%% Set CFGs
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
    latency = cfg.latency;
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

% Set Parameter
if isfield(cfg,'parameter')
    param = cfg.parameter; % What do we want to work on
else
    param = 'avg';
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
Mat = zeros(nsub*ncond,4); % Subjects X Conditons

% Groups
Gr = [];
for g=1:ngr
    tmp = ones(sgr(g)*ncond,1)*g;
    Gr = [Gr; tmp ];
end
Mat(:,2) = Gr; % Group

% Conditions
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
Mat(:,4) =  Subj;

%% Plug in Data
% Look for the right element and put it into the first column of Mat
tic
% Should we average over channels?
if isfield(cfg,'avgoverchan')
    if strcmpi(cfg.avgoverchan,'yes')
        for v = 1:length(varargin)
            varargin{v}.(param) = mean(varargin{v}.(param)(inchan,:));
        end
        nchan = 1; % Update Number of Channels to 1!
        inchan = 1; % Update Channel Input to 1
    end
end

% Should we average over time
if isfield(cfg,'avgovertime')
    if strcmpi(cfg.avgovertime,'yes')
        for v = 1:length(varargin)
            varargin{v}.(param) = mean(varargin{v}.(param)(:,timerange(1):timerange(end)),2);
            varargin{1}.time = 1;
        end
        latency = 1;
    end
end
%% Let's go!
for n = 1:nchan % For Channels
    fprintf('Statistical Stuff at Channel %i of %i \n', n, nchan)
    for t = 1:length(latency);
    fprintf('Statistical Stuff at Timepoint %i of %i \n', t, length(latency))
   
    % Build the interweaved order for the ANOVA
    varind = 0; % Set index to 0

    for g = 1:ngr % Group-loop
        a = 1:sgr(g); % Build a first vector for the index of the first condition
        b = [];       % Build n vectors for the next conditions
        tmp = [];
        for nc = 1:ncond-1
            tmp = (nc*sgr(g)+1):((nc+1)*sgr(g));
            b = [b;tmp];
        end
        c = [a;b]; % concatenate the condition vectors
        tmp = c(:); % flatten
        c = tmp + varind(end); % count up 
        varind = [varind;c]; % stitch together
    end
    varind=varind(2:end);
     
    fi = nearest(varargin{1}.time,latency(t));
    for s = 1:length(varind)
        Mat(s,1) = varargin{varind(s)}.(param)(inchan(n),fi);
    end
    dat = Mat(:,1);

    %% Descriptive Stats
    for g = 1:ngr
        for c = 1:ncond % Freq in this case!
            % Output is ChannelXGroupXFreq
            means(n,t,g,c) = squeeze(mean(dat((Mat(:,3)==c) & (Mat(:,2)==g) ,:)));
            sems(n,t,g,c) = squeeze(std(dat((Mat(:,3)==c) & (Mat(:,2)==g) ,:)))/(sqrt(length(dat((Mat(:,3)==c) & (Mat(:,2)==g) ,:))));
        end
    end

    %% Compute the Anova
    % Blantantly used the existing function but commented the output out
    %  Trujillo-Ortiz
    out = vt_BWAOV2(Mat,alpha); % Hacked Anova Function
    
     %% Put into little boxes
    
    stats.statgroup(n,t,:) = out.F1;
    stats.probgroup(n,t,:) = out.P1;
    stats.dfgroup(:,n) = [out.v1, out.v2];
    stats.statcond(n,t,:) = out.F2;
    stats.probcond(n,t,:) = out.P2;
    stats.dfcond(:,n) = [out.v3, out.v2];
    stats.statint(n,t,:) = out.F3;
    stats.probint(n,t,:) = out.P3;
    stats.dfint(:,n) = [out.v4, out.v5];
    end
end
   
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
            for t = 1:size(stats.statgroup,2)
               % fprintf('Checking Channel %d \n',c)
                if stats.probgroup(c,t,:) < alpha % Group
                    nb =  cellfun(@str2num,neigh(c).neighblabel);
                    % Get Neighbour Indices
                    for cn = 1:length(nb)
                        a(cn) = find(nb(cn) == lab);
                    end
                    % Sort the Indexed Channels
                    nb_alpha = sort(stats.probgroup(a,t)); 
                    if nb_alpha(minnb) < alpha
                        stats.maskgroup(c,t) = 1;
                    else
                        stats.maskgroup(c,t) = 0;
                    end % If Neighbour-Alpha
                end % If Group

                clear a nb nb_alpha

                if stats.probcond(c,t,:) < alpha % Condition
                    nb =  cellfun(@str2num,neigh(c).neighblabel);
                    % Get Neighbour Indices
                    for cn = 1:length(nb)
                        a(cn) = find(nb(cn) == lab);
                    end
                    % Sort the Indexed Channels
                    nb_alpha = sort(stats.probcond(a,t)); 
                    if nb_alpha(minnb) < alpha
                        stats.maskcond(c,t) = 1;
                    else
                        stats.maskcond(c,t) = 0;
                    end % If Neighbour-Alpha
                end % If Cond

                clear a nb nb_alpha

                if stats.probint(c,t,:) < alpha % Interaction
                    nb =  cellfun(@str2num,neigh(c).neighblabel);
                    % Get Neighbour Indices
                    for cn = 1:length(nb)
                        a(cn) = find(nb(cn) == lab);
                    end
                    % Sort the Indexed Channels
                    nb_alpha = sort(stats.probint(a,t)); 
                    if nb_alpha(minnb) < alpha
                        stats.maskint(c,t) = 1;
                    else
                        stats.maskint(c,t) = 0;
                    end % If Neighbour-Alpha
                end % If Interaction

                clear a nb nb_alpha
            end % For Time
        end % For Channel

end % Switch
toc

%% add FT-Stuff
stats.dimord = varargin{1}.dimord;
stats.time = latency;
stats.label = varargin{1}.label(inchan);

