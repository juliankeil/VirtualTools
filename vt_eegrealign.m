function elec_new = vt_eegrealign(cfg,vol,elec_old)

% This function align the standard EEG-electrode positions automatically to
% the BEM-Model. It is based on the nut_eegmralign-function by Dan Wong &
% Sarang Dalal.
% Attention! It needs the Helsinki BEM-Library!
%
% Use as:
% cfg.channel = which channels to include (e.g. 1:126)
% cfg.plot = 'yes' or 'no'
% vol = BEM-Model as obtained by ft_prepare_headmodel
% elec_old = standard electrode positions
%
% Output is a fieldtrip-compatible electrode structure
%
% (c) Julian Keil, 30.07.2014
% Version 1.: First implementation

%%

% Check if PrepareTriangleMesh is in path

if ~exist('PrepareTriangleMesh.m','file')
    errordlg('This function requires the Helsinki BEM library in your Matlab path. It can be downloaded from http://peili.hut.fi/BEM')
    return
end

% Check number of channels
if isfield(cfg,'channel')
    elec_old.chanpos = elec_old.chanpos(cfg.channel,:);
end

% check format of headmodel positions
if isfield(vol.bnd(1), 'pos')
  posfield = 'pos';
elseif isfield(vol.bnd(1), 'pnt')
  posfield = 'pnt';
end

%% Check Units!
if ~strcmpi(vol.unit,elec_old.unit)
    elec_old = ft_convert_units(elec_old, vol.unit); % Make sure its all the same!
end
%% Align using Dan's code
% First we're getting a bit more info on the mesh, i.e. midpoints etc.
mrifid.hsp = PrepareTriangleMesh(vol.bnd(1).(posfield),vol.bnd(1).tri);
  
% Then we'll loop through the channels
for i = 1:length(elec_old.chanpos)
     % First we find the triangle closest to our channel
     senstrimap = dsearchn(mrifid.hsp.mp,elec_old.chanpos(i,:));
     % Then we'll compute the distance to the surface and set the new
     % position
     mrSensorCoord(i,:) = elec_old.chanpos(i,:) - ...
        mrifid.hsp.n(senstrimap,:).*dot(mrifid.hsp.n(senstrimap,:),elec_old.chanpos(i,:)-mrifid.hsp.mp(senstrimap,:))/sum(mrifid.hsp.n(senstrimap,:).^2);
end

%% A plot for checkup
if isfield(cfg,'plot')
    if strcmpi(cfg.plot,'yes')
        figure
        triplot(vol.bnd(3).(posfield), vol.bnd(3).tri,  [], 'faces_skin');
        hold
        triplot(vol.bnd(1).(posfield), vol.bnd(1).tri,  [], 'edges');
        plot3(elec_old.elecpos(:,1),elec_old.elecpos(:,2),elec_old.elecpos(:,3),'r.','MarkerSize',20)
        plot3(mrSensorCoord(:,1),mrSensorCoord(:,2),mrSensorCoord(:,3),'b.','MarkerSize',20)
        camlight left
        view([-90 0]);
    end
end

%% And some cleanup

elec_new = elec_old;
elec_new.chanpos = mrSensorCoord;
elec_new.elecpos = mrSensorCoord;
elec_new.label = elec_old.label(cfg.channel);

