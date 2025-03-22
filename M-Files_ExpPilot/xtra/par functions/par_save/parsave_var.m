% Use this function to save structure data while using parfor
% 
function parsave_var(fname, idx)
  save(fname, 'idx')
end
