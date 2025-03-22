function Surf_fin = smooth_transition(A,B,L_Trans,index_L,index_R)

%%% Input
% Surf1: Main Surface (PIVSurf)
% Surf2: Additional Surface (LFV_Surface_Resized)
% L_Trans = length transition (500)
% index_L = left position of the external matching on Surf1
% index_R = right position of the external matching on Surf1
% L_Trans = 1000;
% index_P = 2000;

w1L = 1:-1/(L_Trans-1):0;
w2L = 0:1/(L_Trans-1):1;
a = A(1:index_L-L_Trans-50);
bL = w1L.*A(index_L-50-L_Trans+1:index_L-50)+w2L.*B(index_L-50-L_Trans+1:index_L-50);
C = B(index_L-50+1:index_R+50);
dR = w1L.*B(index_R+50+1:index_R+50+L_Trans)+w2L.*A(index_R+50+1:index_R+50+L_Trans);
e = A(index_R+50+L_Trans+1:end);
Surf_fin = [a,bL,C,dR,e];
Surf_fin(index_L-100:index_L-1) = smooth(Surf_fin(index_L-100:index_L-1),0.01);
Surf_fin(index_R+1:index_R+100) = smooth(Surf_fin(index_R+1:index_R+100),0.01);