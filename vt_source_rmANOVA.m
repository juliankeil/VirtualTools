function stats = vt_source_rmANOVA(cfg,varargin)
%% ANOVA TEST
% Compute 2-Factor Within Subject Condition and Interactions for all Source
% Points. Use with non-interpolated source analysis files. Make sure that
% all source files have the same number of grid points!
%
% Watch out! You need to apply template grid positions aftewards, as the
% function simply uses the grid positions of the first input!
%
% cfg.channel = Grid indices/names
% cfg.avgoverchan = Average over Grid, or compute voxel-wise ANOVA
% cfg.IV1 = number of Levels in Factor 1
% cfg.IV2 = number of Levels in Factor 2
% cfg.parameter
% cfg.neighbours = neighbourhood structure for spatial correction
% cfg.correctm = 'neighbours' or 'fdr';
% cfg.minnb = minimum number of neighbours, default = 1;
%
% USE AS:
% cfg=[];
% cfg.nIV1 = 2;
% cfg.nIV2 = 2;
% cfg.parameter = 'pow';
% cfg.alpha = .001;
% 
% stats = vt_source_rmANOVA(cfg,fft_vf{:},fft_tf{:},fft_vv{:},fft_tv{:});
%
%
% If you want to use the neighbourhood correction, you have to fake an
% electrode structure and add the "label" subfield of the elec structure to
% the input source files.
%
% elec.elecpos = source.pos;
% elec.chanpos = elec.elecpos;
% for j = 1:length(source.pos);
%    elec.label{j} = sprintf('%d',j);
% end
% 
% cfg=[];
% cfg.elec = elec;
% cfg.method = 'distance';
% cfg.neighbourdist = 2.3;
% 
% neigh = ft_prepare_neighbours(cfg);
%
% (c) Julian Keil 02.06.2015
% Version 1.0.: Implemented Function based on the other ANOVA functions
% (12.12.2014)
% Version 2.: Bugfix in the preparation for the Anova-Matrix. Now supports
% nXm-Anovas (15.06.2017)

%% Set CFGs

nIV1 = cfg.nIV1; % Factors on Level1
nIV2 = cfg.nIV2; % Factors on Level2
nsub = numel(varargin)/(nIV1*nIV2); % How many subjects in total?
subind = 1:nsub;

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

if isfield(cfg,'parameter')
    param = cfg.parameter; % What do we want to work on
else
    param = 'pow';
end

% Set ROI Options
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
    nchan = size(varargin{1}.avg.(param),1);
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
            varargin{v}.avg.(param) = mean(varargin{v}.avg.(param)(inchan,:));
        end
        nchan = 1; % Update Number of Channels to 1!
        inchan = 1; % Update Channel Input to 1
    end
end

%% Let's go!
for n = 1:nchan % For Channels
    fprintf('Statistical Stuff at Source %i of %i \n', n, nchan)
    dat = [];
    tmp = [];

    for c=1:nIV1*nIV2 % Loop Conditions
        for s=1:nsub % loop Subjects varagin
            sid = sub2ind([nsub, nIV1*nIV2],s,c); % Subject Identifier, Subjects are order 1...n;1...n.. in varargin 
            tmp_g(sid,1) = varargin{sid}.avg.(param)(inchan(n));  % take the relevant element from the input
        end % for Conditions  
    end % for Subjects

    Mat(:,1) = tmp_g;
    AnovaMat(:,:) = Mat;

    %% Descriptive Stats
    
    for g = 1:nIV1
        for c = 1:nIV2 % 
            % Output is ChannelXIV1XIV2
            tmp = nanmean(AnovaMat((AnovaMat(:,3)==c) & (AnovaMat(:,2)==g),:));
            means(n,g,c) = tmp(1); 
            tmp_z = nanstd(AnovaMat((AnovaMat(:,3)==c) & (AnovaMat(:,2)==g) ,:));
            tmp_n = sqrt(length(AnovaMat(AnovaMat(:,3)==c & AnovaMat(:,2)==g)));
            sems(n,g,c) = tmp_z(1)/tmp_n(1);
        end
    end
        
    stats.means = means;
    stats.sems = sems;
    %% Compute the Anova
    % Blantantly used the existing function but commented the output out
    %  Trujillo-Ortiz
        out = vt_RMAOV2(squeeze(AnovaMat(:,:)),alpha); % Hacked Anova Function

    %% Put into little boxes

    stats.statIV1(n,:) = out.F1;
    stats.probIV1(n,:) = out.P1;
    stats.dfIV1(:,n) = [out.v1, out.v2];
    stats.statIV2(n,:) = out.F2;
    stats.probIV2(n,:) = out.P2;
    stats.dfIV2(:,n) = [out.v3, out.v2];
    stats.statint(n,:) = out.F3;
    stats.probint(n,:) = out.P3;
    stats.dfint(:,n) = [out.v5, out.v6];
    
end

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

        [stats.probIV1_fdr(:,:), stats.maskIV1(:,:)] = vt_fdr(stats.probIV1(:,:), alpha);
        [stats.probIV2_fdr(:,:), stats.maskIV2(:,:)] = vt_fdr(stats.probIV2(:,:), alpha);
        [stats.probint_fdr(:,:), stats.maskint(:,:)] = vt_fdr(stats.probint(:,:), alpha);

    case{'neighbours'}
        fprintf('Using Neighbours to correct for for multiple comparisons \n')

            stats.maskIV1(:,:) = zeros(size(stats.statIV1(:,:)));
            stats.maskIV2(:,:) = zeros(size(stats.statIV2(:,:)));
            stats.maskint(:,:) = zeros(size(stats.statint(:,:)));
            % Get Input Labels
            if iscell(varargin{1}.label)
                lab = varargin{1}.label;
            else
                lab = cellfun(@str2num,varargin{1}.label); 
            end
            
            for c = 1:length(stats.statIV1(:,:))
               % fprintf('Checking Channel %d \n',c)
                if stats.probIV1(c,:) < alpha % Group
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
                    nb_alpha = sort(stats.probIV1(a,:)); 
                    if nb_alpha(minnb) < alpha
                        stats.maskIV1(c,:) = 1;
                    else
                        stats.maskIV1(c,:) = 0;
                    end % If Neighbour-Alpha
                end % If Group

                clear nb nb_alpha

                if stats.probIV2(c,:) < alpha % Condition
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
                    nb_alpha = sort(stats.probIV2(a,:)); 
                    if nb_alpha(minnb) < alpha
                        stats.maskIV2(c,:) = 1;
                    else
                        stats.maskIV2(c,:) = 0;
                    end % If Neighbour-Alpha
                end % If Cond

                clear nb nb_alpha

                if stats.probint(c,:) < alpha % Interaction
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
                    nb_alpha = sort(stats.probint(a,:)); 
                    if nb_alpha(minnb) < alpha
                        stats.maskint(c,:) = 1;
                    else
                        stats.maskint(c,:) = 0;
                    end % If Neighbour-Alpha
                end % If Interaction

                clear nb nb_alpha

            end % For Channel
end % Switch
toc

%% add FT-Stuff
stats.inside = varargin{1}.inside;
stats.outside = varargin{1}.outside;
stats.pos = varargin{1}.pos;
stats.dim = varargin{1}.dim;

