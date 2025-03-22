%%% Video writer

v = VideoWriter([TestName,'_solo_onde_da_vento_StereoPIV_dati_base.avi']); % apri il file video
open(v)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    my_vertices=[[X(end,1), X(1,:), X(end,end), X(end,1)]; [Y(1,1), Lev, Y(1,1), Y(1,1)]]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


f=figure(1);
f.Position = [100 100 1240 840];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disegna le velocità medie di fase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
qq=quiver(X(1:end,1:end),Y(1:end,1:end),Ufase(1:end,1:end)/20,Vfase(1:end,1:end)/20,0,'linewidth',2);
qq.Color='[0.8 0.8 0.8]';

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disegna le velocità istantanee al netto delle velocità medie di fase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
qq=quiver(X(1:end,1:end),Y(1:end,1:end),(Ucorr(1:end,1:end)-Ufase(1:end,1:end))/20,(Vcorr(1:end,1:end)-Vfase(1:end,1:end))/20,0,'linewidth',2);
qq.Color='r';

qq1=quiver(xmax-60,ymax-15,10,0,0,'linewidth',2);
qq1.Color='r';
text(xmax-70,ymax-10,'20 cm/s', 'FontSize', 24)

rr=plot(X(1,:),Lev,'linewidth',2);
rr.Color='b';



set(findall(gcf,'-property','FontSize'),'FontSize',24)
ax = gca;
ax.LineWidth = 3;


frame = getframe(gcf);
writeVideo(v,frame)
hold off

close(v)
