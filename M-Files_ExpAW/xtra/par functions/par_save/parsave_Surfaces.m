% Use this function to save structure data while using parfor
% 
function parsave_Surfaces(fname, PixRes, PIVRes)
  save(fname, 'PixRes', 'PIVRes')
end
