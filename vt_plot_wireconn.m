function jk_plot_wireconn(cfg,stats)

% Plots a wire frame of the brain along with the color coded sources of
% interest
%
% cfg.paramter = plot parameter
% cfg.dummy = dummy structure with positions of grid points
% cfg.vol = volume of brain -> BEM Model
% cfg.freq = Frequency
% cfg.source = Label of Source
% cfg.sink = Label of Sink
% cfg.sourceplot = option to plot Ortho-Plot as well. 'yes' or 'no'
%                  If sourceplot = 'yes', you also have to specify an MRI
% cfg.mri = MNI-MRI (e.g. T1.nii from SPM)
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
% Version 1.3.: November 2013: 
%                               Added Support for multiple Sinks

freq=nearest(stats.freq,cfg.freq);
vol=cfg.vol;
dummy = cfg.dummy;
parameter = cfg.parameter;
indata = stats.(parameter);
sourcelab = str2num(cfg.source{1});

% Loop through the sinks
for s = 1:length(cfg.sink)
    sinklab{s} = str2num(cfg.sink{s});
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

% Then plot
% First Check for FT-Version
% Plot the Wiremesh
if exist('vt_triplot') == 0 % Check which version of triplot is available
    triplot(vol.bnd(3).pos, vol.bnd(3).tri,  [], 'edges');
else
    vt_triplot(vol.bnd(3).pos, vol.bnd(3).tri,  [], 'edges');
end

hold; % Hold the Wiremesh

% Plot the Grid
plot3(dummy.pos(dummy.inside,1),dummy.pos(dummy.inside,2),dummy.pos(dummy.inside,3),'b*','LineWidth',1,'Color',[.5 .5 .5]);
plot3(dummy.pos(sourcelab,1),dummy.pos(sourcelab,2),dummy.pos(sourcelab,3),'g*','LineWidth',5);

% And now plot the connectivity lines
for s=1:length(sinklab)
    plot3(dummy.pos(sinklab{s},1),dummy.pos(sinklab{s},2),dummy.pos(sinklab{s},3),'r*','LineWidth',5);
    line([dummy.pos(sourcelab,1) dummy.pos(sinklab{s},1)],...
    [dummy.pos(sourcelab,2) dummy.pos(sinklab{s},2)],...
    [dummy.pos(sourcelab,3) dummy.pos(sinklab{s},3)],'LineWidth',5,'Color',hot2(end-s,:));
end

% Plot the Colorbar
colormap(hot3);
hcb=colorbar;
set(hcb,'YTick',[0.00000001 .5 1],'YTickLabel',[mat_s(1) mat_s(ceil(length(mat_s)/2)) mat_s(end)]);

%% Should there be an ortho-plot as well?

if strcmpi(cfg.sourceplot,'yes');
    
    % prepare source-plot
    mri = cfg.mri;
    dummy.avg.pow2(dummy.inside) = 0;
    dummy.avg.pow2(dummy.outside) = NaN;
    dummy.avg.pow2(sourcelab) = mean(mat_s);
    for s=1:length(sinklab)
        dummy.avg.pow2(sinklab{s}) = mat_s(end-s);
    end
    
    tmpcfg = [];
    tmpcfg.parameter = 'avg.pow2';
    tmpcfg.downsample = 1;
    
    dummy_i = ft_sourceinterpolate(tmpcfg,dummy,mri);
    
    % And Plot the Ortho
    tmpcfg = [];
    tmpcfg.method='ortho';
    tmpcfg.funparameter = 'avg.pow2';
    tmpcfg.maskparameter = tmpcfg.funparameter;
    
    ft_sourceplot(tmpcfg,dummy_i);
    
end