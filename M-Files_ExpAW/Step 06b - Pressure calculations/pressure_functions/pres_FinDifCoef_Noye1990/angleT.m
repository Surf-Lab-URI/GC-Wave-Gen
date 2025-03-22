function a = angleT(x,y)

if x>0
    a = atan(y/x);
else
    a = pi+atan(y/x);
end