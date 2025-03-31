clear
clc
close all
%%
load('GravityWaveSurfacesSample.mat')
%%
figure(3)
plot(minPts(:,1), minPts(:,2),'.r','MarkerSize',15)
%%
clear pts pt lastPts nPts
pt = minPts(1,:);
pts(1,:) = pt;
nPts = 1;
lastPts = [];
lines = [];

for i = 2:length(minPts)
    if pt(2) == minPts(i,2)
        pts(nPts + 1,:) = minPts(i,:);
        nPts = nPts+1;
    else
        if length(lastPts) == 0
            lastPts = pts;
            nPts = 0;
            pt = minPts(i,:);
            pts = [];
            pts(1,:) = pt;
        else
            nLines = length(lines);
            % nStepLines = 0;
            steplines = [];
            stepLines = table(steplines);
            for k = 1:length(lastPts)
                for l = 1:length(pts)
                    stepline = {[lastPts(k,:);pts(l,:)]};
                    stepLines = [stepLines;table(stepline,'VariableNames','steplines')];
                end
            end
        end
    end

    pt = minPts(i,:);
end
