function vt_plot_weathermap(cfg,stats)

% Plots Weathermap for Stats
% Compute proper FT-Stats first
%
% inputs
% cfg.parameter = parameter to plot
% cfg.order = 'yes' or 'no' - should the plot be ordered for regions
% cfg.computemask = 'yes' or 'no' - compute mask for input values on the
% fly
% if yes, then set:
% cfg.numconseq = how many timepoints should be significant (default = 10)
% cfg.alpha = alpha level (default = .05)
% cfg.probparameter = parameter containing P-Values
% cfg.maskparameter = Previously computed Mask (e.g. Neighbour or
% FDR-Correction)
%
% if no, then you can set:
% cfg.maskparameter = maskparameter
%
% (c) Julian Keil 2014 - Based on Daniel Senkowski's lodi weathermap plot
% Version 1.0.: 5.3.2014
% Version 1.1.: 5.3.2014: Stupid bugfixes for masking and ordering
% Version 1.2.: 30.05.2016: Added integration with precomputed masks
% Version 1.3.: 20.06.2016: Added switch from logical to double for the
% precomputed masks
% Version 1.3.1.: 22.06.2016: Added proper color bars for one-sided
% distributions

%% Set cfgs

parameter = cfg.parameter;
x_axis = stats.time;

%% Change Mask from Logical to Double if necessary
if islogical(stats.maskIV1)
    stats.maskIV1 = double(stats.maskIV1);
    stats.maskIV2 = double(stats.maskIV2);
    stats.maskint = double(stats.maskint);
end

%% Reorder Channels?

if isfield(cfg,'order')
    if strcmpi(cfg.order,'yes')
        new_order = [25 21 16 17 22 26 29 30 27 23 18 13 12 8 9 14 19 24 28 31 32 83 77 15 10 6 5 3 65 66 68 72 78 88 93 84 73 96 67 70 74 79 89 85 80 75 71 76 81 86 90 94 95 91 87 82 92 96 ,...
                    60 56 20 51 50 55 59 63 62 58 54 49 54 46 11 7 41 44 48 53 43 40 37 4 2 1 34 36 39 47 57 61 64 52 42 38 35 33 97 98 101 106 111 121 116 110 105 100 99 102 103 104 109 115 120 126 125 119 114 108 107 112 113 118 124 123 117 122];
        % Assign stat values
        val = stats.(parameter)(new_order,:);
        
        if isfield(cfg,'probparameter')
          stats.(cfg.probparameter) = stats.(cfg.probparameter)(new_order,:);
        end
        
        if isfield(cfg,'maskparameter')
          stats.(cfg.maskparameter) = stats.(cfg.maskparameter)(new_order,:);
        end
    else
        val = stats.(parameter);
    end % if order = yes
    
else
    val = stats.(parameter);
end % if order

%% Set the Zeros in the Binary Mask to NaNs to avoid artifacts with Prob-Parameter
if isfield(cfg,'maskparameter')
    for x = 1:length(stats.(cfg.maskparameter)(:))
        if stats.(cfg.maskparameter)(x) == 0;
            stats.(cfg.maskparameter)(x) = NaN;
        end
    end
    if isfield(cfg,'probparameter')
       stats.(cfg.probparameter) = stats.(cfg.probparameter).* stats.(cfg.maskparameter);
    end
end
%% Mask stat values for consequtive timepoints

if isfield(cfg,'computemask') % Check for Masking
    if strcmpi(cfg.computemask,'yes') % Should we compute a mask?
        if isfield(cfg,'numconseq')
          numconseq = cfg.numconseq;
        else
          numconseq = 10;
        end % if numcons

        if isfield(cfg,'alpha')
            alphalevel = cfg.alpha;
        else
            alphalevel = .05;
        end % if alpha

        if isfield(cfg,'probparameter')
          prob = stats.(cfg.probparameter);
        else
          display('Please set P-Value-Parameter')
          return
        end % if prob

        % Actually Compute the Mask
        T = zeros(size(val)); % New Variable to put masked values in
            for i=1:size(val,1) % For each channel
                for k=1:size(val,2) % For each Timepoint
                    if (k - numconseq) > 0 % Check if there is enough time left
                        if prob(i,k-[0:numconseq]) < alphalevel % Check if Value is smaller than alpha
                            T(i,k-[0:numconseq])=val(i,k-[0:numconseq]); % Set to Value
                        else
                            T(i,k)=0; % Set to 0
                        end % If alpha
                    end % If Time
                end % For Time
            end % For Chan
            
    elseif isfield(cfg,'maskparameter') % Should we use a precomputed mask?
            T = val.* stats.(cfg.maskparameter);
    else % Let's not mask at all
            T = val;
    end % If Compute
      
end % If Mask
    
%% Open up the figure
figure
        
p = imagesc(x_axis,[1:size(T,1)], T);
axis tight

if min(T(:)) == 0
    caxis([0 max(abs(T(:)))]);
    colorbar;

    % Build a custom colorbar
        rval = ones(100,1);
        gval = (0:.01:.99)';
        bval = (0:.01:.99)';

        zero = [1 1 1];

        cmap = flipud([rval gval bval; zero]);
          
elseif max(T(:)) == 0;
    caxis([-max(abs(T(:))) 0]);
    colorbar;

    % Build a custom colorbar
        rval = ones(100,1);
        gval = (0:.01:.99)';
        bval = (0:.01:.99)';

        zero = [1 1 1];

        cmap = flipud([zero; flipud(bval) flipud(gval) flipud(rval)]); 
else  
    caxis([-max(abs(T(:))) max(abs(T(:)))]);
    colorbar;

    % Build a custom colorbar
        rval = ones(100,1);
        gval = (0:.01:.99)';
        bval = (0:.01:.99)';

        zero = [1 1 1];

        cmap = flipud([rval gval bval; zero; flipud(bval) flipud(gval) flipud(rval)]);

end
colormap(cmap);
%

if isfield(cfg,'order')
    if strcmpi(cfg.order,'yes')
        set(gca,'ytick',[3 20 46 75 101 119]);
        set(gca, 'YTickLabel', 'fronto-polar|frontal|fronto-central|central|parieto-occipital|occipital');
    end
end
        
set(gca,'FontSize',10)          
          
end
