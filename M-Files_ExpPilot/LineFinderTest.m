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
lines = table(lines);

figure(4)
plot(minPts(1,1),minPts(1,2),'.r','MarkerSize',15)
hold on

for i = 2:length(minPts)
    plot(minPts(i,1),minPts(i,2),'.r','MarkerSize',15)
    hold on
    if pt(2) == minPts(i,2)
        pts(nPts + 1,:) = minPts(i,:);
        nPts = nPts+1;
    else
        if length(lastPts) == 0
            lastPts = pts;
            nPts = 1;
            pt = minPts(i,:);
            pts = [];
            pts(1,:) = pt;
        else
            nLines = height(lines);
            nStepLines = 0;
            steplines = [];
            stepLines = table(steplines);
            for k = 1:size(lastPts,1)
                for l = 1:size(pts,1)
                    steplines = {[lastPts(k,:);pts(l,:)]};
                    stepLines = [stepLines;table(steplines)];
                    nStepLines = nStepLines+1;
                    clear ll
                    ll = cell2mat(stepLines.steplines(nStepLines));
                    plot(ll(:,1), ll(:,2),'-b')
                end
            end
 
            if nLines > 0
                
            end

        end
    end

    pt = minPts(i,:);
end
