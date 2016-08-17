function stat = vt_compare_corr(cfg,varargin);

% Compare Different Correlation Values
%
% This Function is an implementation of Eq. 6.95 and 6.96 from Bortz
%
% Input:
% cfg.alpha = alpha level, default 0.05
% cfg.rvals = r values to compare
% cfg.sgr = group size, vector with group sizes belonging to rvals 
% cfg.correctm = correction method, right now only 'neighbours'
% cfg.neighbours = neighbourhood-structure to use with correction
% cfg.latency = time of interest
% cfg.channel = channels of interest
% Output:
% * V = Chi2-Distributed test statistic
% * crtival = critical Chi2-Value
% * mask = 0 if V below critval, 1 if above
%
% Ver 1.: Julian Keil, 08.01.2015

%% Set defaults
tic
stat = [];

if isfield(cfg,'alpha')
    p = 1-alpha;
else
    p = .95;
end
% Number of Groups
ngr = length(cfg.sgr);
df = ngr-1;

% Select Parameter

if isfield(cfg,'param')
    param = cfg.param;
else
    param = 'cor';
end

% Set Correction 
if isfield(cfg,'correctm') % Check for cfg
    correctm = cfg.correctm; % set Method
    if strcmpi(correctm,'none')
        correctm = 'none';
        fprintf('Not using correction for multiple comparisons \n')
    elseif strcmpi(correctm,'neighbours') && isfield(cfg,'neighbours')
        neigh = cfg.neighbours; % Use Predefined Neighbours
        fprintf('Using Neighbour Correction \n')
    else
        fprintf('Please provide Neighbourhood Structure \n');
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
    %ntime = length(cfg.latency); % How many Time Steps
    %tstep = varargin{1}.time(2)-varargin{1}.time(1);
    timerange = [nearest(varargin{1}.time,cfg.latency(1)):1:nearest(varargin{1}.time,cfg.latency(end))]; 
    ntime = length(timerange);
    end
else       
    ntime = length(varargin{1}.time); % How many Time Steps
    timerange = 1:ntime;    
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
            tmp =  cellfun(@num2str,num2cell(cfg.channel),'UniformOutput',false);
            for c = 1:length(tmp)
                inchan(c) = find(strcmpi(tmp{c},varargin{1}.label)==1);
            end
        end
        
    end
else
    nchan = size(varargin{1}.(param),1);
    inchan = 1:nchan;
end
%% Loop Channels and TimePoints
disp('Comparing correlations...')
for c = 1:size(varargin{1}.(param),1)
    for t = 1:size(varargin{1}.(param),2)
        
        for i = 1:length(varargin)
            dat(i) = varargin{i}.(param)(c,t);
        end
        % Zscore r-avlues

        Zj = 0.5.*(log(1+dat)-log(1-dat));

        %% First Compute U

        U = sum((cfg.sgr-ngr).*Zj)/sum(cfg.sgr-ngr); % Eq. 6.96

        %% Second, compute V

        V = sum((cfg.sgr-ngr).*(Zj-U).^2);

        %% Finally, compare V to Chi2-Table

        critval = chi2inv(p,df); % Look up Chi2-Table

%         if V >= critval
%             stat.mask(c,t) = 1;
%         else
%             stat.mask(c,t) = 0;
%         end

        stat.V(c,t) = V;
        stat.critval(c,t) = critval;
    end % Time
end % Channel

%% Correction
switch correctm
    case{'none'}
        fprintf('Not performing correction for multiple comparisons \n Mask is uncorrected! \n')
            stat.mask = stat.V>=stat.critval;
case{'neighbours'}
        fprintf('Using Neighbours to correct for for multiple comparisons \n')
        for f = 1:ntime
            stat.mask(:,f,:) = zeros(size(stat.V(:,f,:)));

            % Get Input Labels
            if iscell(varargin{1}.label)
                lab = varargin{1}.label;
            else
                lab = cellfun(@str2num,varargin{1}.label); 
            end
            
            for c = 1:length(stat.V(:,f,:))
               % fprintf('Checking Channel %d \n',c)
                if stat.V(c,f,:) >= critval % Group
                    if iscell(neigh(c).neighblabel)
                        nb = neigh(c).neighblabel;
                        % Get Neighbour Indices
                        for cn = 1:length(nb)
                            a(cn) = find(strcmpi(nb(cn),lab));
                        end
                    else
                        nb =  cellfun(@str2num,neigh(c).neighblabel);
                        % Get Neighbour Indices
                        for cn = 1:length(nb)
                            a(cn) = find(nb(cn) == lab);
                        end
                    end
                    
                    % Sort the Indexed Channels
                    nb_alpha = sort(stat.V(a,f,:)); 
                    if nb_alpha(minnb) >= critval
                        stat.mask(c,f,:) = 1;
                    else
                        stat.mask(c,f,:) = 0;
                    end % If Neighbour-Alpha
                end % If Group

                clear nb nb_alpha

            end % For Channel
        end % For Freq
end % End Switch

stat.time = varargin{1}.time;

toc
