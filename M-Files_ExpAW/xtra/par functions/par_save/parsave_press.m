% Use this function to save structure data while using parfor
% 
function parsave_press(fname, press, Lapl)
  save(fname, 'press', 'Lapl')
end