clear all
close all

fiche=fopen('E:\data\RampeDeVent\LC.txt','r') ;
LC1=fscanf(fiche,'%g %g',[2 inf]);
LC1=LC1';
% LC11=LC1(5001:end,1)-10;
LC11=LC1(15501:end,1)-31;

data=LC1;
calib_Pa = [0.0039 0.0531 0.1047 0.1539 0.2046 0.2542 0.3049 0.3549 0.4042 0.4534 0.5035]*248.84; %in pascal
calib_VDC = [0.0884 0.5759 1.0952 1.5882 2.0989 2.5985 3.1077 3.6072 4.0983 4.5847 5.0815]; %in volts
poly_pit = polyfit(calib_VDC, calib_Pa, 1);
rho = 1.2;
uPit = (2*(polyval(poly_pit,data)/rho));
L = length(data);
u = real(sqrt(uPit));
LC1 = u;
fclose(fiche);

fiche=fopen('E:\data\RampeDeVent\LC2.txt','r') ;
LC2=fscanf(fiche,'%g %g',[2 inf]);
LC2=LC2';
data=LC2;
calib_Pa = [0.0039 0.0531 0.1047 0.1539 0.2046 0.2542 0.3049 0.3549 0.4042 0.4534 0.5035]*248.84; %in pascal
calib_VDC = [0.0884 0.5759 1.0952 1.5882 2.0989 2.5985 3.1077 3.6072 4.0983 4.5847 5.0815]; %in volts
poly_pit = polyfit(calib_VDC, calib_Pa, 1);
rho = 1.2;
uPit = (2*(polyval(poly_pit,data)/rho));
L = length(data);
u = real(sqrt(uPit));
LC2 = u;
fclose(fiche);

fiche=fopen('E:\data\RampeDeVent\LC3.txt','r') ;
LC3=fscanf(fiche,'%g %g',[2 inf]);
LC3=LC3';
LC22=LC2(5001:end,1)-10;
data=LC3;
calib_Pa = [0.0039 0.0531 0.1047 0.1539 0.2046 0.2542 0.3049 0.3549 0.4042 0.4534 0.5035]*248.84; %in pascal
calib_VDC = [0.0884 0.5759 1.0952 1.5882 2.0989 2.5985 3.1077 3.6072 4.0983 4.5847 5.0815]; %in volts
poly_pit = polyfit(calib_VDC, calib_Pa, 1);
rho = 1.2;
uPit = (2*(polyval(poly_pit,data)/rho));
L = length(data);
u = real(sqrt(uPit));
LC3 = u;
fclose(fiche);

fiche=fopen('E:\data\RampeDeVent\LC4.txt','r') ;
LC4=fscanf(fiche,'%g %g',[2 inf]);
LC4=LC4';
data=LC4;
calib_Pa = [0.0039 0.0531 0.1047 0.1539 0.2046 0.2542 0.3049 0.3549 0.4042 0.4534 0.5035]*248.84; %in pascal
calib_VDC = [0.0884 0.5759 1.0952 1.5882 2.0989 2.5985 3.1077 3.6072 4.0983 4.5847 5.0815]; %in volts
poly_pit = polyfit(calib_VDC, calib_Pa, 1);
rho = 1.2;
uPit = (2*(polyval(poly_pit,data)/rho));
L = length(data);
u = real(sqrt(uPit));
LC4 = u;
fclose(fiche);

% inter1 = [1:40000];
% inter2 = inter1+5000;
inter1 = [1:29500];
inter2 = inter1+15500;
figure, plot(LC11(inter1,1),LC1(inter2,2))
hold on,  plot(LC11(inter1,1),LC2(inter2,2),'r'), plot(LC11(inter1,1),LC3(inter2,2),'k'), plot(LC11(inter1,1),LC4(inter2,2),'c')

