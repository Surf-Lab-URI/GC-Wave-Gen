function [surf_s] = FiltSurf(surf, tol)
x = 1:length(surf);
f = spaps(x,surf,tol);

surf_s = fnval(f,x);
end