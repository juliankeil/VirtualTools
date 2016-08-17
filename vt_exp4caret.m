function vt_exp4caret(cfg,stats)

% Saves the Grid-Structure Source datat to an .nii file that can then be
% plotted using Caret 
% see here for a tutorial: 
% http://prefrontal.org/blog/2009/04/using-caret-for-fmri-visualization/
%
% cfg.paramter = plot parameter
% cfg.dummy = dummy structure with positions of grid points
% cfg.vol = volume of brain -> BEM Model
% cfg.freq = Frequency
% cfg.toi = Time of Interest
% cfg.mri = MNI-MRI (e.g. T1.nii from SPM)
% cfg.outname = name to save to
%
% (c) Julian Keil, 2014
% Version 1.0.: 02.04.2014
% Version 1.1.: 25.06.2014: Added Freq Range Option

%% Set Basics
vol=cfg.vol;
dummy = cfg.dummy;
parameter = cfg.parameter;
indata = stats.(parameter);
outname = cfg.outname;

if isfield(cfg,'freq')
    freq=nearest(stats.freq,cfg.freq);
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

%% 
mat = indata(:,freq);
ia = cellfun(@str2num,stats.label);

%% prepare source-data
mri = cfg.mri;
if length(mat) < length(dummy.inside)
    for c = 1:length(ia)
        [a(c), b(c)] = find(dummy.inside == ia(c));
    end
    dummy.avg.pow2(dummy.inside(b)) = mat;
else
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
    
%% interpolate Source Data
tmpcfg = [];
tmpcfg.parameter = 'avg.pow2';
tmpcfg.downsample = 1;

dummy_i = ft_sourceinterpolate(tmpcfg,dummy,mri);

if isfield(dummy_i,'time')
    dummy_i = rmfield(dummy_i,'time');
end

%% Export as .nii
cfg=[];
cfg.parameter = 'avg.pow2';
cfg.filename = outname;
cfg.datatype = 'float';

ft_volumewrite(cfg,dummy_i);