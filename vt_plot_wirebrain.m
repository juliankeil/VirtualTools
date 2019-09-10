function vt_plot_wirebrain(cfg,stats)

% Plots a wire frame of the brain along with the color coded sources of
% interest
%
% cfg.paramter = plot parameter
% cfg.dummy = dummy structure with positions of grid points
% cfg.vol = volume of brain -> BEM Model
% cfg.freq = Frequency
% cfg.toi = Time of Interest
% cfg.sourceplot = option to plot Ortho-Plot as well. 'yes' or 'no'
%                  If sourceplot = 'yes', you also have to specify an MRI
% cfg.mri = MNI-MRI (e.g. T1.nii from SPM)
% cfg.islocation = 'yes' or 'no': Were the virtual channels built on
% isolated locations?
% cfg.atlas = path to atlas
%
% (c) Julian Keil 2013
%
% Version 1.0.: June 2013
% Version 1.1.: October 2013: 
%                               Added switch for ft_triplot vs. triplot
%                               Changed "dummy.stats.pos" to "dummy.pos"
%                               Changed the color value assignment
% Version 1.2.: October 2013:
%                               Added Ortho-Plot Option
% Version 1.3.: October 2013:   
%                               Added Time-Option
% Version 1.4.: June 2014:      Added Frequency-Band Option
% Version 1.5.: July 2015:      Added Atlas Option


vol=cfg.vol;
dummy = cfg.dummy;
parameter = cfg.parameter;
indata = stats.(parameter);

% check format of headmodel positions
if isfield(vol.bnd(1), 'pos')
  posfield = 'pos';
elseif isfield(vol.bnd(1), 'pnt')
  posfield = 'pnt';
end

if isfield(cfg,'freq')
    freq=nearest(stats.freq,cfg.freq(1));
    if length(cfg.freq) > 1
        f1 = nearest(stats.freq,cfg.freq(1));
        f2 = nearest(stats.freq,cfg.freq(end));
        indata = mean(indata(:,f1:f2),2);
        freq=1;
    end
elseif isfield(cfg,'toi')
    t1 = nearest(stats.time,cfg.toi(1));
    t2 = nearest(stats.time,cfg.toi(2));
    toi = [t1 t2];
    indata = mean(indata(:,toi),2);
    freq = 1; % Dummy Value
end

% Set location option
if ~isfield(cfg, 'islocation') || strcmp(cfg.islocation,'no')
  islocation = 0;
elseif strcmp(cfg.islocation,'yes')
  islocation = 1;
else
  error('cfg.islocation should be set to ''yes'' or ''no''.')
end

% Set plot option
if ~isfield(cfg, 'sourceplot') || strcmp(cfg.sourceplot,'no')
  sourceplot = 0;
elseif strcmp(cfg.sourceplot,'yes')
  sourceplot = 1;
else
  error('cfg.sourceplot should be set to ''yes'' or ''no''.')
end

% Set Atlas
if isfield(cfg,'atlas');
    atlas = cfg.atlas;
else
    atlas = [];
end

%% Make a Color Map
hot2=hot(length(indata(:,freq)));

%% Plot Color Coded Stats Values
% First sort Color Map by Stat-VAlues
mat = indata(:,freq);
[mat_s, im] = sort(mat); % get the data index
hot3 = hot2; % Save colors for later

% set zero values to zero ;-)
for m =1:length(mat_s)
    if mat_s(m) == 0 || mat_s(m) == 1
        hot2(m,:) = [.5 .5 .5]; % Zero values are grey
    end
end

% Then get the Source IDs

if islocation
    ia = 1:length(stats.label);
    ia = ia(im);
else
    ia = cellfun(@str2num,stats.label);
    ia = ia(im); % sort the IDs by data index
end

% Then plot
% First Check for FT-Version
% The Plot the Wiremesh
if exist('vt_triplot') == 0 % Check which version of triplot is available
    triplot(vol.bnd(3).(posfield), vol.bnd(3).tri,  [], 'edges');
else
    vt_triplot(vol.bnd(3).(posfield), vol.bnd(3).tri,  [], 'edges');
end

hold; % Hold the Wiremesh

% And now plot the indata values to each point
for c=1:length(mat_s)
    plot3(dummy.pos(ia(c),1),dummy.pos(ia(c),2),dummy.pos(ia(c),3),'*','Color',hot2(c,:),'LineWidth',5);
end

% Plot the Colorbar
colormap(hot3);
hcb=colorbar;
set(hcb,'YTick',[0.00000001 .5 1],'YTickLabel',[mat_s(1) mat_s(ceil(length(mat_s)/2)) mat_s(end)]);

%% Should there be an ortho-plot as well?

if sourceplot
    
    % prepare source-plot
    mri = cfg.mri;
    if length(mat) < length(dummy.inside)
        for c = 1:length(ia)
            [a(c), b(c)] = find(dummy.inside == ia(c));
        end
        dummy.avg.pow2(dummy.inside(b)) = mat;
    else
        dummy.avg.pow2 = dummy.avg.pow;
        dummy.avg.pow2(dummy.inside) = mat;
    end
        
    dummy.avg.pow2(dummy.outside) = NaN;
    
    % Flip Dimensions if neccessary (possibly when using location grid)
    if size(dummy.avg.pow) ~= size(dummy.avg.pow2) 
        dummy.avg.pow2 = dummy.avg.pow2';
    end
    
    % Add dimensions if neccessary
    if ~isfield(dummy,'dim')
        dim1 =  length(min(dummy.pos(:,1)):max(dummy.pos(:,1)));
        dim2 =  length(min(dummy.pos(:,2)):max(dummy.pos(:,2)));
        dim3 =  length(min(dummy.pos(:,3)):max(dummy.pos(:,3)));
        
        dummy.dim = [dim1 dim2 dim3];
    end
    
    tmpcfg = [];
    tmpcfg.parameter = 'pow2';
    tmpcfg.downsample = 1;
    
    dummy_i = ft_sourceinterpolate(tmpcfg,dummy,mri);
    
    if isfield(dummy_i,'time')
        dummy_i = rmfield(dummy_i,'time');
    end
    
    % And Plot the Ortho
    tmpcfg = [];
    tmpcfg.method='ortho';
    %tmpcfg.surffile = 'surface_l4_both.mat';
    tmpcfg.funparameter = 'pow2';
    tmpcfg.maskparameter = tmpcfg.funparameter;
    tmpcfg.atlas = atlas;
    
    ft_sourceplot(tmpcfg,dummy_i);
    
end
