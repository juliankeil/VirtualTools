function stats = vt_time_rmANOVA_tmp(cfg,varargin)
%% ANOVA TEST
% Compute 2-Factor Within Subject Condition and Interactions for all EEG-Sensors
%
% cfg.channel = channel indices/names
% cfg.avgoverchan = Average over Channels, or compute channel-wise ANOVA
% cfg.avgovertime = Average over Time
% cfg.latency = Latency to compute the ANOVA for;
% cfg.IV1 = number of Levels in Factor 1
% cfg.IV2 = number of Levels in Factor 2
% cfg.parameter = Input data field
% cfg.correctm = 'none', 'cluster', 'neighbours' or 'fdr';
% cfg.neighbours = neighbourhood structure for spatial correction
% cfg.minnb = minimum number of neighbours, default = 1;
% cfg.numrandomization = number of randomizations; Set to at least 2.

%
% USE AS:
% cfg=[];
% cfg.latency = [.1 .15];
% cfg.avgoverfreq = 'yes';
% cfg.nIV1 = 2;
% cfg.nIV2 = 2;
% cfg.parameter = 'avg';
% cfg.alpha = .001;
% 
% stats = vt_time_rmANOVA(cfg,fft_vf{:},fft_tf{:},fft_vv{:},fft_tv{:});
%
% (c) Julian Keil 16.04.2014
% Version 1.0.: Implemented Function
% Version 1.1.: Adapted for Timelock Data
% Version 1.2.: Corrected Error in Averaging
% Version 2.0.: Updated Channel Selection Option (04.11.2014)
% Version 2.1.: Updated Channel Selection for Neighbour-Correction
% (12.12.2014)
% Version 3.: MAjor Buugfix in Neighbour-Correction (04.11.2015)
% Version 3.1.: More Bugfixes in Correction thanks to Mathis (30.05.2016)
% Version 4.: Bugfix in the preparation for the Anova-Matrix. Now supports
% nXm-Anovas (15.06.2017)
% Version 4.1.: Removed unneccessary time-loops.
% Version 5.: Added cluster correction option. This option will randomly
% shuffle the condition labels and comute the size of the largest resulting
% cluster. All clusters in the original data below the alpha-level
% perceptile will be set to 0. (jk, 29.08.2019)

%% Set CFGs

nIV1 = cfg.nIV1; % Factors on Level1
nIV2 = cfg.nIV2; % Factors on Level2
nsub = numel(varargin)/(nIV1*nIV2); % How many subjects in total?
%subind = 1:nsub;

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
    elseif strcmpi(correctm,'neighbours') || strcmpi(correctm,'cluster') && isfield(cfg,'neighbours')
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
if isfield(cfg,'latency')
    if isfield(cfg,'avgovertime')
        if strcmpi(cfg.avgovertime,'yes')
            ntime = 1;
            timerange = [nearest(varargin{1}.time,cfg.latency(1)) nearest(varargin{1}.time,cfg.latency(end))];
        end
    else
    %ntime = length(cfg.latency); % How many Time Steps
    %tstep = varargin{1}.time(2)-varargin{1}.time(1);
    timerange = nearest(varargin{1}.time,cfg.latency(1)):1:nearest(varargin{1}.time,cfg.latency(end)); 
    ntime = length(timerange);
    end
else       
    ntime = length(varargin{1}.time); % How many Time Steps
    timerange = 1:ntime;    
end

if isfield(cfg,'parameter')
    param = cfg.parameter; % What do we want to work on
else
    param = 'avg';
end

% Set Channel Options
if isfield(cfg,'channel')
    if strcmpi(cfg.channel,'all')
        nchan = size(varargin{1}.(param),1);
        inchan = 1:nchan;
    else
        nchan = length(cfg.channel);
        if iscell(cfg.channel)
            inchan = zeros(1,length(cfg.channel));
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
for s=1:nIV1*nIV2
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
            varargin{v}.(param) = mean(varargin{v}.(param)(inchan,:));
        end
        nchan = 1; % Update Number of Channels to 1!
        inchan = 1; % Update Channel Input to 1
    end
end

% Should we average over frequencies
if isfield(cfg,'avgovertime')
    if strcmpi(cfg.avgovertime,'yes')
        for v = 1:length(varargin)
            varargin{v}.(param) = mean(varargin{v}.(param)(:,timerange(1):timerange(end)),2);
            varargin{1}.time = 1;

        end
        %latency = 1;
    end
end

%% Let's go!
% If Cluster correction, then find the clusters belonging together
if strcmpi(correctm,'cluster')
    nrand = cfg.numrandomization;
    
    % Build SpatDimNeighStruct for findclusters
    dummymat = zeros(size(neigh,2),size(neigh,2));
    for c = 1:size(neigh,2)
        dummymat(c,c) = 1;
        for n = 1:size(neigh(c).neighblabel,1)
            for ic = 1:size(neigh,2)
                if strcmp(neigh(ic).label,neigh(c).neighblabel(n))
                    dummymat(c,ic) = 1;
                end % if compare labels
            end % For All Channels
        end % for Neighbours
    end % For Channels
else % Else just do one loop
    nrand = 1;
end % If cluster

% Pre-allocate variables
clusterIV1 = zeros(1,nrand);
clusterIV2 = zeros(1,nrand);
clusterint = zeros(1,nrand);

means = NaN(nchan,ntime,nIV1,nIV2);
sems = means;

for r = 1:nrand
    randvec = randperm(length(Mat(:,1)));
    fprintf('Statistical Stuff: Permutation %i \n', r)
    for n = 1:nchan % For Channels
        %dat = [];
        %tmp = [];
        if ntime == 1
            timerange = 1;
        end

        for f = 1:ntime % loop frequencies
            fi = timerange(f);
            for c=1:nIV1*nIV2 % Loop Conditions
                for s=1:nsub % loop Subjects varagin
                    sid = sub2ind([nsub, nIV1*nIV2],s,c); % Subject Identifier, Subjects are order 1...n;1...n.. in varargin 
                    tmp_g(sid,1) = varargin{sid}.(param)(inchan(n),fi);  % take the relevant element from the input
                end % for Conditions  
            end % for Subjects

            Mat(:,1) = tmp_g;

        %% Descriptive Stats
            for g = 1:nIV1
                for c = 1:nIV2 % 
                    % Output is ChannelXFreqXIV1XIV2
                    tmp = mean(Mat((Mat(:,3)==c) & (Mat(:,2)==g),:));
                    means(n,f,g,c) = tmp(1); 
                    tmp_z = std(Mat((Mat(:,3)==c) & (Mat(:,2)==g) ,:));
                    tmp_n = sqrt(length(Mat(Mat(:,3)==c & Mat(:,2)==g)));
                    sems(n,f,g,c) = tmp_z(1)/tmp_n(1);
                end
            end

        %% Compute the Anova
        % Blantantly used the existing function but commented the output out
        %  Trujillo-Ortiz
            if r == 1
                out = vt_RMAOV2(Mat,alpha); % Hacked Anova Function

                % Put into little boxes
                stats.statIV1(n,f,:) = out.F1;
                stats.probIV1(n,f,:) = out.P1;
                stats.dfIV1(:,n,f) = [out.v1, out.v2];
                stats.statIV2(n,f,:) = out.F2;
                stats.probIV2(n,f,:) = out.P2;
                stats.dfIV2(:,n,f) = [out.v3, out.v2];
                stats.statint(n,f,:) = out.F3;
                stats.probint(n,f,:) = out.P3;
                stats.dfint(:,n,f) = [out.v5, out.v6];
                
            else
                tmpval = Mat(randvec,1); % apply the same random permutation to all values in one shuffle
                tmpMat = Mat;
                tmpMat(:,1) = tmpval;
                tmpout = vt_RMAOV2(tmpMat,alpha); % Hacked Anova Function
                
                tmpstat.P1(n,f,:) = tmpout.P1;
                tmpstat.P2(n,f,:) = tmpout.P2;
                tmpstat.P3(n,f,:) = tmpout.P3;
                
            end % r            
        end % time
    end % channel
    
    % If we're in the randomization loop, find the largest clusters
    if r > 1
        % Find the Clusters per Permutation.
        maskIV1 = tmpstat.P1<alpha;
        maskIV2 = tmpstat.P2<alpha;
        maskint = tmpstat.P3<alpha;
        
        clear cluster1 cluster2 cluster3
        [cluster1, num] = vt_findcluster(maskIV1,dummymat,minnb);
        [cluster2, num] = vt_findcluster(maskIV2,dummymat,minnb);
        [cluster3, num] = vt_findcluster(maskint,dummymat,minnb);

        if any(cluster1(:))
            clusterIV1(r) = sum(cluster1(:) == 1);
        else
            clusterIV1(r) = 0;
        end

        if any(cluster2(:))
            clusterIV2(r) = sum(cluster2(:) == 1);
        else
            clusterIV2(r) = 0;
        end

        if any(cluster3(:))
            clusterint(r) = sum(cluster3(:) == 1);
        else
            clusterint(r) = 0;
        end % if cluster
    end % if r   
end % For Nrand

stats.means = means;
stats.sems = sems;
stats.dimord = 'chan_time';

%% Correction
switch correctm
    case{'none'}
        fprintf('Not performing correction for multiple comparisons \n Mask is uncorrected! \n')
            stats.maskIV1 = stats.probIV1<alpha;
            stats.maskIV2 = stats.probIV2<alpha;
            stats.maskint = stats.probint<alpha;
            
    case{'cluster'} % Compare empirical clusters to shuffled clusters
        fprintf('Using clusters based on shuffled data to estimate cluster size \n Mask is uncorrected, Cluster-Masks contain Cluster infos \n')
            
            % Mask P-Values
            stats.maskIV1 = stats.probIV1<alpha;
            stats.maskIV2 = stats.probIV2<alpha;
            stats.maskint = stats.probint<alpha;
            
            % Find Clusters in P-Values
            [cluster1, num1] = vt_findcluster(stats.maskIV1,dummymat,minnb);
            [cluster2, num2] = vt_findcluster(stats.maskIV2,dummymat,minnb);
            [cluster3, num3] = vt_findcluster(stats.maskint,dummymat,minnb);
            
            % Find Alpha-Percentile in Shuffled Clusters
            y1 = quantile(clusterIV1,1-alpha);
            y2 = quantile(clusterIV2,1-alpha);
            yint = quantile(clusterint,1-alpha);
            
            % Set all Empirical Clusters below Alpha-Percentile to 0
            if num1 > 0
                for cl = 1:num1
                    if sum(cluster1(:) == cl) <= y1
                        cluster1(cluster1 == cl) = 0; % Remove the too-small cluster
                    end
                end
            else
                fprintf('No Clusters for IV1 found \n')
            end
            stats.clustermaskIV1 = cluster1;
            
            if num2 > 0
                for cl = 1:num2
                    if sum(cluster2(:) == cl) <= y2
                        cluster2(cluster2 == cl) = 0; % Remove the too-small cluster
                    end
                end
            else
                fprintf('No Clusters for IV2 found \n')
            end
            stats.clustermaskIV2 = cluster2;
            
            if num3 > 0
                for cl = 1:num3
                    if sum(cluster3(:) == cl) <= yint
                        cluster3(cluster3 == cl) = 0; % Remove the too-small cluster
                    end
                end
            else
                fprintf('No Clusters for int found \n')
            end
            stats.clustermaskint = cluster3;
            
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
            for t = 1:ntime
                    
                if stats.probIV1(c,t) < alpha % is smaller than alpha?
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
                    nb_alpha = sort(stats.probIV1(a,t)); 
                    if nb_alpha(minnb) < alpha % is the nth neighbour also smaller alpha?
                        stats.maskIV1(c,t) = 1;
                    else
                        stats.maskIV1(c,t) = 0;
                    end % If Neighbour-Alpha
                end % If Group

                clear a nb nb_alpha

                if stats.probIV2(c,t) < alpha % Condition
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
                    nb_alpha = sort(stats.probIV2(a,t)); 
                    if nb_alpha(minnb) < alpha
                        stats.maskIV2(c,t) = 1;
                    else
                        stats.maskIV2(c,t) = 0;
                    end % If Neighbour-Alpha
                end % If Cond

                clear a nb nb_alpha

                if stats.probint(c,t) < alpha % Interaction
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
                    nb_alpha = sort(stats.probint(a,t)); 
                    if nb_alpha(minnb) < alpha
                        stats.maskint(c,t) = 1;
                    else
                        stats.maskint(c,t) = 0;
                    end % If Neighbour-Alpha
                end % If Interaction

                clear a nb nb_alpha

            end % For Time
        end % For Chan
end
toc

%% add FT-Stuff
stats.time = varargin{1}.time(timerange);
stats.label = varargin{1}.label(inchan);

