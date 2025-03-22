% Use this function to save structure data while using parfor
% 
function parsave_CartVel(fname, Cartesian, PixRes, PIVRes, CST)
  save(fname, 'Cartesian', 'PixRes', 'PIVRes', 'CST')
end
