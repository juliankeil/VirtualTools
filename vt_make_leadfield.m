function [lf] = vt_make_leadfield(cfg)

% Builds a leadfield based on a given MRI and different input parameters.
% Options are:
%   - Create a new BEM-Model with specified resolution based on a MRI
%   - Create a leadfield with a specified resolution based on a Head-Model
%   - Create a leadfield for a subset of grid points based on a ROI

% (c) Julian Keil 2013
% Ver 0.1.: 10.09.2013

%%

