function data=trialchooser_new(data, nep)
%data=trialchooser(data, nep)
%chooses randomnly n-epochs from preprocessing structure.

vec=1:length(data.trial);
vec2=vec(randperm(length(data.trial)));

if length(data.trial)<nep
    error('There are not enough trials.')
end

ind=vec2(1:nep);

data.trial=data.trial(ind);
data.time=data.time(ind);

% if isfield(data.cfg, 'trl')
%    data.cfg.trl=data.cfg.trl(ind,:);
% end

if isfield(data, 'sampleinfo')
    data.sampleinfo = data.sampleinfo(ind,:);
end

if isfield(data, 'trialinfo')
    data.trialinfo = data.trialinfo(ind,:);
end