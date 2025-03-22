% function [L] = par_load(LaplName)
%   
%   A = load(LaplName);
%   AA = fields(A);
%   L = eval(['A.' AA{1} ';']);
% end

function [varargout] = par_load(structName)
  
  A = load(structName);
  AA = fields(A);
  varargout{1} = eval(['A.' AA{1} ';']);
  if numel(AA)>1
      for i = 2:numel(AA)
          varargout{i} = eval(['A.' AA{i} ';']);
      end
  end
