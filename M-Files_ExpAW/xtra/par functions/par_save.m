% Use this function to save structure data while using parfor 
% 
function par_save(varargin)
savefile = varargin{1}; % first input argument
for i=2:nargin
    savevar.(inputname(i)) = varargin{i}; % other input arguments
end
save(savefile,'-struct','savevar')

%%% ALERT: gives transparency error when used in parfor. 
%%% To save the pressure, use par_save_press
%%% To save the Laplacian, use par_save_lapl