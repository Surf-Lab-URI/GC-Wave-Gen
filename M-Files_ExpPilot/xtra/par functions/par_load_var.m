% function [L] = par_load(LaplName)
%   
%   A = load(LaplName);
%   AA = fields(A);
%   L = eval(['A.' AA{1} ';']);
% end

% Load variables from a single structure

function [varargout] = par_load_var(structName,varargin)
  
    varargout = cell(length(varargin),1);
    for i = 1:length(varargin)
        varargout{i} = load(structName,varargin{i});
    end