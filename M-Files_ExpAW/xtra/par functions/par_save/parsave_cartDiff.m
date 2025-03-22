% Use this function to save structure data while using parfor
% 
function parsave_cartDiff(fname, cartDiff, PixRes, PIVRes)
  save(fname, 'cartDiff', 'PixRes', 'PIVRes')
end
