function lf_roi = vt_make_roifield(cfg)

% Builds a Leadfield for ROIs defined in MNI-Space
% Input needs a MNI-normalized MRI, an atlas or a mask 
% 
% cfg.mri = mni-normalized MRI
% cfg.roi = Name of Roi or MNI location in case of spherical ROI
% 
% You can chose the ROI either based on an atlas OR a logical mask OR on a
% single location and sphere radius
% cfg.atlas = where to find the atlas
% OR
% cfg.mask = whole mri structure from http://neuro.imm.dtu.dk/services/jerne/ninf/voi.html
% If you use a mask from this site, read the .hdr and .img-file with 
% mask = ft_read_mri('nameofthefile.hdr')
% OR
% cfg.sphere = radius of ROI, unit is the same as in cfg.lf
% OR
% cfg.box = Nx3 Definition of Box around location
% 
% cfg.vol = mni-normalized Headmodel
% cfg.lf = leadfield based on vol
% cfg.hemisphere = left right or both
%
% WATCH OUT! You need SPM(8) to correctly read in the .hdr and .img-files
%
%
% Use as:
%cfg=[];
%cfg.mri = mri; % Standard MRI
%cfg.vol = vol; % Individual Headmodel
%cfg.standardvol = standardvol; % Standard Headmodel for plotting
%cfg.lf = lf_126; % Individual Sourcemodel in Standard Space!
%cfg.roi = 'Superior Temporal Gyrus';
%cfg.atlas = '/home/keil/data/programs/Matlab/fieldtrip-20130618/template/atlas/afni/TTatlas+tlrc.BRIK';
%cfg.hemisphere = 'left';
%lf_new = vt_make_roifield(cfg);

%
% (c) Julian Keil 2013
% Version 1.0.: 12.09.2013
% Version 1.1.: 22.10.2013: Added Mask-Option (Martin und Julian)
% Version 1.2.: 23.10.2013: Bug in hemisphere selection fixed (Martin und Julian)
% Version 1.2.1: 21.11.2013: Fixed the same bug as above again.
% Version 1.3.: 21.11.2013: Added Sphere Option
% Version 1.4.: 29.09.2014: Added Box Option
% Version 1.5.: 18.03.2015: Changed the use of triplot to ft_plot_vol
% Version 1.6.: 28.09.2016: Fixed an error in the Pythagoras Formula
% Version 1.7.: 31.05.2023: Updated to new use of ft_prepare_sourcemodel
% and added wait bar

%% Set cfgs
% Load MRI
mri = cfg.mri;

        % Downsample MRI
        xfg = [];
        xfg.downsample = 2;

        mri = ft_volumedownsample(xfg,mri);
        
% Load Standard Volume and convert to mm
standardvol = cfg.standardvol;
standardvol = ft_convert_units(standardvol,'mm');
        
lf = cfg.lf;
vol = cfg.vol;
if isfield(cfg,'roi') && isfield(cfg,'atlas') % ROI based on atlas
    ROIdef = cfg.roi;
    roitype = 1;
    atlas = cfg.atlas;
    
elseif isfield(cfg,'roi') && isfield(cfg,'sphere')  % Spherical Roi with Center on 'roi' 
                                                % and radius as 'sphere'
    ROIdef = cfg.roi;
    ROIrad = cfg.sphere;
    roitype = 3;
    
elseif isfield(cfg,'roi') && isfield(cfg,'box')  % Box Roi with Center on 'roi' 
                                                % and radius as 'sphere'
    ROIdef = cfg.roi;
    ROIrad = cfg.box;
    roitype = 4;

elseif isfield(cfg,'mask')
    tmpcfg=[];
    tmpcfg.parameter = 'anatomy';
    tmp = ft_sourceinterpolate(tmpcfg,cfg.mask,mri);
    
    ROIdef = tmp;
    roitype = 2;
end

if isfield(cfg,'hemisphere')
    hem = cfg.hemisphere;
else
    hem = 'both';
end


%% Set the start and End of the x,y,z axis
vox1 = mri.transform*[mri.dim,1]'; % start
vox1 = vox1(1:3)/10; % convert to mm
vox2 = mri.transform*[1,1,1,1]';   %end
vox2 = vox2(1:3)/10; % convert to mm

%% Choose the Minima
vox_min = min([vox1,vox2]')';
vox_max = max([vox1,vox2]')';

%% Set the spacing
vox_delta = abs(mri.transform*[1,1,1,1]'-mri.transform*[2,2,2,1]'); 
vox_delta = vox_delta(1:3)/10; % convert to mm
RES = vox_delta(1);

% Now we have the grid definition for the MNI Brain -> Super dense!

%% Build a grid on the MNI voxels
% this dummy gradiometer definition is required for prepare_dipole_grid and
% headmodelplot. We don't actually need a leadfield, just the grid
% positions!
template_grad     = [];
template_grad.pnt = [];
template_grad.ori = [];
template_grad.tra = [];
template_grad.label = {};
template_grad.unit = 'cm';

% construct the dipole grid in the template brain coordinates
% the source units are in cm
% the negative inwardshift means an outward shift of the brain surface for
% inside/outside detection
cfg = [];
%cfg.inwardshift = -0.5;
cfg.method = 'basedongrid';
cfg.xgrid = vox_min(1):RES:vox_max(1); 
cfg.ygrid = vox_min(2):RES:vox_max(2); 
cfg.zgrid = vox_min(3):RES:vox_max(3);
cfg.tight = 'no';
cfg.headmodel = vol;
cfg.grad = template_grad;

template_grid = ft_prepare_sourcemodel(cfg);

%% Convert to mm if lf is in mm
if strcmpi(lf.unit,'mm')
    template_grid = ft_convert_units(template_grid,'mm');
end

%% Get the ROI based on atlas
switch roitype
    case(1)
        cfg=[];
        cfg.atlas = atlas;
        %cfg.inputcoord = 'mni';
        cfg.roi = ROIdef;

        ROI = ft_volumelookup(cfg,mri);
    
    case(2)
        ROI = ROIdef.anatomy;  
        
    case(3)
        cfg=[];
        %cfg.inputcoord = 'mni';
        cfg.roi = ROIdef;
        cfg.sphere = ROIrad; % Radius of Sphere

        ROI = ft_volumelookup(cfg,mri);
        
    case(4)
        cfg=[];
        %cfg.inputcoord = 'mni';
        cfg.roi = ROIdef;
        cfg.box = ROIrad; % Radius of Sphere

        ROI = ft_volumelookup(cfg,mri);
    
end
%% Now the funky part: Use 3D-Pythagoras s=sqrt((X2-X1)²+(Y2-Y1)²+(Z2-Z1)²)
% to compute the distance between the template_grid and the pre-computed leadfield
f = waitbar(0,'Calculating distance...');
l = length(lf.inside);
data = zeros(l,1);

for i=1:length(lf.inside)
    %fprintf('Getting the position of Point %i of %i\n',i,length(lf.inside))
    waitbar(i/length(lf.inside),f,sprintf('Getting the position of Point %i of %i',i,l));
    vox = lf.inside(i);
    dist = sqrt(sum((ones(size(template_grid.pos,1),1)*lf.pos(vox,:) - template_grid.pos ).^2,2));
    [dummy,index(i)] = min(dist);
    data(i) = ROI(index(i));
end
close(f)

%% Some Clean-Up
% divide all points in data by maximum within data -> Normalize to [0 1]
% data=data/max(data);
% % set all data < 0.01 to zero
% data(find(data<0.01)) = 0;

%%
% find all indices of data points > 0
index = find(data>0.1);

mask = zeros(size(data));
mask(index)=1;

%% find position of index points within the leadfield

switch hem
    case{'both'} 
        fprintf('ROI is on both hemispheres')
        temp = lf.pos(lf.inside(index),:);
        index_ROI = lf.inside(index);
        mask_ROI = zeros(size(data));
        mask_ROI(index_ROI) = 1;
    case{'left'}
        fprintf('ROI is on left hemisphere')
        temp = lf.pos(lf.inside(index),:);
        temp_i =  index(temp(:,1)<0);
        index_ROI = lf.inside(temp_i);
        mask_ROI = zeros(size(data));
        mask_ROI(index_ROI) = 1;
    case{'right'}
        fprintf('ROI is on right hemisphere')
        temp = lf.pos(lf.inside(index),:);
        temp_i =  index(temp(:,1)>0);
        index_ROI = lf.inside(temp_i);
        mask_ROI = zeros(size(data));
        mask_ROI(index_ROI) = 1;
end

%%
if isfield(cfg,'plot')
    if strcmpi(cfg.plot,'yes')
        figure;
        ft_plot_headmodel(standardvol, 'edgecolor', 'none'); alpha 0.3;
        hold
        plot3(lf.pos(index_ROI,1),lf.pos(index_ROI,2),lf.pos(index_ROI,3),'r.','MarkerSize',20)
        camlight left
        view([-90 0]);
    end
end
 
 %% Output
 
 lf_roi = lf;
 lf_roi.inside = index_ROI;
 lf_roi.outside = setdiff(1:length(lf.pos),lf_roi.inside);
