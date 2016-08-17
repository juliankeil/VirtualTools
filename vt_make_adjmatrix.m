function [outdata] = vt_make_adjmatrix(cfg,indata)
%
% Build Adjacency Matrix for Graphtheoretical Measures
% i.e. mask out all connection below threshold
%
% Input:
% 
% cfg.method = 'proportion', 'fixed' or 'relative'
%
% WATCH OUT! If you want to compare different conditions/groups/time
% points, you should use the same threshold in all adjacency matrices!
%
% cfg.threshold = either a fixed value or x times std
% Default relative is above 1.5 std
% Proportional threshold can be between 1 (all connections preserved) and 0
% (no connections preserved)
%
% cfg.metric = connectivity metric
% Can be 'plv' 'coh' or 'psi' or you can add your method of choice below
% line 38
%
% (c) Julian Keil, 2013
% Version 1: 27.11.2013
% Version 2: 05.12.2013 - Added Proportional Thresholding
% Version 3: 14.02.2014 - Fixed Thresholding for Bidirectional data (Upper
% and lower cutoff)
% Version 4: 27.06.2014 - Added Z-Score Threshold Option: Statistically
% significant values are above 1.645 (one-sided test on absolute values
% p<.05)
% Version 5: 12.01.2015 - Added Frequency Specfic  Fixed Value Option
% Version 6: 25.07.2016 - Replaced NaNs with 0s
%% Set cfgs

if strcmpi(cfg.method, 'fixed') && isfield(cfg,'threshold')
    thresh = cfg.threshold;
elseif strcmpi(cfg.method, 'fixed') && ~isfield(cfg,'threshold')
    fprint('Please enter threshold or use relative method')
    return
elseif strcmpi(cfg.method, 'relative') && isfield(cfg,'threshold')
    thresh = cfg.threshold;
elseif strcmpi(cfg.method, 'fixed') && ~isfield(cfg,'threshold')
    thresh = 1.5;
elseif strcmpi(cfg.method,'proportion') && isfield(cfg,'threshold')
    thresh = cfg.threshold;
elseif strcmpi(cfg.method,'proportion') && ~isfield(cfg,'threshold')
    fprint('Please enter threshold or use relative method')
    return
elseif strcmpi(cfg.method,'specific') && isfield(cfg,'threshold');
    thresh = cfg.threshold;
elseif strcmpi(cfg.method,'specific') && ~isfield(cfg,'threshold')
    fprint('Please enter threshold or use relative method')
    return
end

if ~isfield(cfg,'metric')
    fprintf('Please indicate the connectivity metric used')
    return
elseif strcmpi(cfg.metric,'plv');
    tmpdata = indata.plvspctrm;
elseif strcmpi(cfg.metric,'coh');
    tmpdata = indata.cohspctrm;
elseif strcmpi(cfg.metric,'psi');
    tmpdata = indata.psispctrm;
elseif strcmpi(cfg.metric,'wpli_debiased');
    tmpdata = indata.wpli_debiasedspctrm;
else
    fprintf('dude, can you please add your metric to this list?\n')
    return
end

%% Compute Mask
% First replace NaNs with zeros
tmpdata(isnan(tmpdata)) = 0;

% Then compute Masks
fprintf('Masking Connectivity...\n')
switch(cfg.method)
    case{'fixed'} % Absolute Threshold
        for f = 1:length(indata.freq)
            for n = 1:length(indata.label)
                tmp = tmpdata(n,:,f); % Select Connections between nodes for freq
                if isfield(cfg,'zscore')
                    if strcmpi(cfg.zscore,'yes')
                        tmpind = find(abs(zscore(tmp)) >= thresh); % Find all elements larger than fixed threshold
                    
                    elseif strcmpi(cfg.zscore,'no')
                        tmpind = find(abs(tmp) >= thresh); % Find all elements larger than fixed threshold
                        
                    end
                        
                else
                        tmpind = find(abs(tmp) >= thresh); % Find all elements larger than fixed threshold
                end
                tmpzeros = zeros(size(tmp));
                tmpzeros(tmpind) = 1;
                adjmatspctrm(n,:,f) = tmpzeros;
            end % nodes
        end % freq
        
    case{'relative'} % Relative Threshold
        for f = 1:length(indata.freq)
            for n = 1:length(indata.label)
                tmp = tmpdata(n,:,f); % Select Connections between nodes for freq
                tmpstd = std(tmp);
                tmpind = find(abs(tmp) >= thresh*tmpstd); % Find all elements larger that 1.5 std
                tmpzeros = zeros(size(tmp));
                tmpzeros(tmpind) = 1;
                adjmatspctrm(n,:,f) = tmpzeros;
            end % nodes
        end % freq
        
    case{'proportion'} % Proportional Threshold. 
        for f = 1:length(indata.freq)
            for n = 1:length(indata.label)
                tmp = tmpdata(n,:,f); % Select Connections between nodes for freq
                [tmps tmpind] = sort(tmp,2,'descend'); % Sort tmps
                tmpzeros = zeros(size(tmp)); % Make a zero-matrix
                
                n_cut = floor(length(tmp)*thresh/2); % How man elements should we trow out?
                tmpu = tmps(1:n_cut); % Throw away these elements, preserve the rest
                tmpl = tmps(length(tmps)-(n_cut-1):end); % Throw away these elements, preserve the rest
                tmpzeros(tmpind(1:n_cut))=tmpu; % And add back into the zero-amtrix
                tmpzeros(tmpind(length(tmps)-(n_cut-1):end))=tmpl; % And add back into the zero-amtrix
                
                adjmatspctrm(n,:,f) = double(tmpzeros~=0);
            end % nodes
        end % Freq
        
    case{'specific'} % Frequency Specific Thresholds
        if length(indata.freq) == length(thresh)
            for f = 1:length(indata.freq)
                for n = 1:length(indata.label)
                    tmp = tmpdata(n,:,f); % Select Connections between nodes for freq
                    
                    % special chase: if thresh is 0, take elements NOT 0
                    if thresh == 0
                        tmpind = find(abs(tmp) > thresh(f)); % Find all elements larger that 1.5 std
                    else
                        tmpind = find(abs(tmp) >= thresh(f)); % Find all elements larger that 1.5 std
                    end
                    
                    if isfield(cfg,'zscore')
                        if strcmpi(cfg.zscore,'yes')
                            tmpind = find(abs(zscore(tmp)) >= thresh(f)); % Find all elements larger than fixed threshold

                        elseif strcmpi(cfg.zscore,'no')
                            tmpind = find(abs(tmp) >= thresh(f)); % Find all elements larger than fixed threshold
                        end
                    else
                        tmpind = find(abs(tmp) >= thresh); % Find all elements larger than fixed threshold
                    end
                    
                    tmpzeros = zeros(size(tmp));
                    tmpzeros(tmpind) = 1;
                    adjmatspctrm(n,:,f) = tmpzeros;
                end % nodes
            end % freq
        else
            disp('Please provide a threshold for each frequency')
        end % If
        
end % switch

%% Clean up

outdata = indata;
outdata.adjmatspctrm = adjmatspctrm;
outdata.dimord = 'chan_freq';

    