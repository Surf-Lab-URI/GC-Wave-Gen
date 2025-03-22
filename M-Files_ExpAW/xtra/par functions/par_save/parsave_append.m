% Use this function to save structure data while using parfor
% 
function parsave_append(fname, idx)
  save(fname, 'idx', '-append')
end
