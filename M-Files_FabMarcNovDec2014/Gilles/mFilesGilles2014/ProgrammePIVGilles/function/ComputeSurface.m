function [xx, yy, sp] = ComputeSurface(I1, firstFrame, x0, y0, sp)
%perform texture analysis in computing the surface of intersection

I = abs(I1 - medfilt2(I1, [3, 3], 'symmetric'));
nhoodSize = 21;
nhood = true(nhoodSize);
pI = entropyfilt(I, nhood);
entropyI = 1 - (pI - min(min(pI)))/(max(max(pI)) - min(min(pI)));
pI = rangefilt(I, nhood);
rangeI = 1 - (pI - min(min(pI)))/(max(max(pI)) - min(min(pI)));
pI = stdfilt(I, nhood);
stdI = (pI - min(min(pI)))/(max(max(pI)) - min(min(pI)));
pI = colfilt(I, [nhoodSize, nhoodSize], 'sliding', @mean);
meanI = (pI - min(min(pI)))/(max(max(pI)) - min(min(pI)));

A = [entropyI(:), meanI(:), stdI(:), rangeI(:)];
phi = A - repmat(mean(A), size(A, 1), 1);
[u, d, v] = svd(phi'*phi);
p = (u'*A');
pI = reshape(p(1, :), size(I, 1), size(I, 2));
linComb = (pI - min(min(pI)))/(max(max(pI)) - min(min(pI)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SNAKE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%knots used for the spline
knots = 10;

% initialize with hand for the first frame
if (firstFrame == true)
   figure, colormap(gray(256)), imagesc(linComb); hold on; axis image;
   pos = ginput;
%     load firstframe.mat;
    %plot(pos(:, 1), pos(:, 2), 'k*');
    [a, b] = krige(pos(:, 1)', pos(:, 2)', 1, size(I, 2));

    % using a least sqaures spline
    spl2 = spap2(knots, 4, a', b');
    x0 = 1:size(I, 2);
    y0 = spval(spl2, x0);

    %creating the image model energy based on likelihood
    a = linComb;
    [xx, yy] = meshgrid(1:size(I, 2), 1:size(I, 1));
    t = yy < repmat(y0, size(yy, 1), 1);
    a(t) = NaN;
    %estimating the mixture of distributions
    [f, xi] = ksdensity(a(:), 'kernel', 'epanechnikov');
    %using spline to interpolate the pdf
    sp = spapi(5, xi, f);

    close all;
end

ML = reshape(spval(sp, linComb(:)), size(linComb));
MLdensity = colfilt(ML, [nhoodSize, nhoodSize], 'sliding', @mean);

% pos now holds the sampled points
subsample = 64;
snakeIx = x0(1:subsample:end);
snakeIy = y0(1:subsample:end);

%window to perform particle filtering
probWy = 96; probWx = 16; likelihoodWin = 16;
%increase the hypothesis space for particle filter
iter = 10;

for loopcount = 1:2
    %number of snake points
    numP = length(snakeIx);
    state = cell(1, numP);

    % window for the particle filtering for first point
    window = MLdensity(floor(snakeIy(1) - probWy:snakeIy(1) + probWy), ...
        floor(snakeIx(1):snakeIx(1) + probWx));
    tState = ParticleFilter(window, likelihoodWin/2, iter);
    deltax = tState(:, 1) - (size(window, 2) + 1)/2;
    deltay = tState(:, 2) - (size(window, 1) + 1)/2;
    state{1, 1} = [deltax, deltay, tState(:, 3)];

    for n = 2:numP-1
        % window for the particle filtering
        window = MLdensity(floor(snakeIy(n) - probWy:snakeIy(n) + probWy), ...
            floor(snakeIx(n) - probWx:snakeIx(n) + probWx));
        tState = ParticleFilter(window, likelihoodWin, iter);
        deltax = tState(:, 1) - (size(window, 2) + 1)/2;
        deltay = tState(:, 2) - (size(window, 1) + 1)/2;
        state{1, n} = [deltax, deltay, tState(:, 3)];
    end

    % window for the particle filtering for last point
    window = MLdensity(floor(snakeIy(numP) - probWy:snakeIy(numP) + probWy), ...
        floor(snakeIx(numP) - probWx:snakeIx(numP)));
    tState = ParticleFilter(window, likelihoodWin/2, iter);
    deltax = tState(:, 1) - (size(window, 2) + 1)/2;
    deltay = tState(:, 2) - (size(window, 1) + 1)/2;
    state{1, numP} = [deltax, deltay, tState(:, 3)];

    % getting the maximum possible elements
    maxH = 1;
    for n = 1:numP
        if length(state{1, n}(:, 1)) > maxH
            maxH = length(state{1, n}(:, 1));
        end
    end

    posx = zeros(maxH, numP);
    posy = zeros(maxH, numP);
    energy = zeros(maxH, numP);
    size_a = [];
    for n = 1:numP
        a = state{1, n}; size_a = [size_a, size(a, 1)];
        posx(1:size(a, 1), n) = a(:, 1) + snakeIx(n);
        posy(1:size(a, 1), n) = a(:, 2) + snakeIy(n);
        energy(1:size(a, 1), n) = a(:, 3);
    end

    normE =  - (energy ./ repmat(max(energy), size(energy, 1), 1));

    %find the positions where the likelihood was confusing and adjacent points
    conft = find(size_a > 1);
    conf = setdiff(union([conft-1, conft, conft + 1],[]), [0, 1, numP]);

    alpha = 1.5; beta = 0.75; kappa = 0.25;
    %using the subsample as the average distance between points
    d_avg = subsample;

    %holds the curv positions before DP
    curv = cell(maxH, numP);
    
    for snakePt = 1:numP
        %handling the first point
        if snakePt == 1
            for k = 1:size_a(1)
                v_i = [posx(k, 1), posy(k, 1)];
                Eext = normE(k, 1);
                %holds the energy terms from the two sides
                tE = [];
                for right = 1:size_a(2)
                    v_i_r = [posx(right, 2), posy(right, 2)];

                    %using the distance from the point to d_avg
                    Edist = abs(d_avg - norm(v_i_r - v_i))/d_avg;
                    %angle between the two vectors
                    tE = [tE; right, alpha*Edist + kappa*Eext, right];
                end
                minpos = find(min(tE(:, 2)) == tE(:, 2));
                curv{k, 1} = tE(minpos(1), :);
            end
        elseif snakePt == numP
            for k = 1:size_a(numP)
                v_i = [posx(k, numP), posy(k, numP)];
                Eext = normE(k, numP);
                %holds the energy terms from the two sides
                tE = [];
                for left = 1:size_a(numP - 1)
                    v_i_l = [posx(left, numP - 1), posy(left, numP - 1)];
                    %using the distance from the point to d_avg
                    Edist = abs(d_avg - norm(v_i - v_i_l))/d_avg;
                    %angle between the two vectors
                    tE = [tE; left, curv{left, numP - 1}(1, 2) + alpha*Edist + kappa*Eext, 1];
                end
                minpos = find(min(tE(:, 2)) == tE(:, 2));
                curv{k, numP} = tE(minpos(1), :);
            end
        else
            %if the current snake point belongs to the confusion set
            num = (snakePt == conf);
            % there exists atleast one confusion state
            if sum(num) >= 1
                for k = 1:size_a(conf(num))
                    v_i = [posx(k, conf(num)), posy(k, conf(num))];
                    Eext = normE(k, conf(num));
                    %holds the energy terms from the two sides
                    tE = [];
                    for left = 1:size_a(conf(num) - 1)
                        v_i_l = [posx(left, conf(num) - 1), posy(left, conf(num) - 1)];
                        for right = 1:size_a(conf(num) + 1)
                            v_i_r = [posx(right, conf(num) + 1), posy(right, conf(num) + 1)];

                            %using the distance from the point to d_avg
                            Edist = abs(d_avg - (norm(v_i_r - v_i) + norm(v_i - v_i_l))/2)/d_avg;
                            % akgul's paper and mapping to [-1 0]
                            Esmo = ((1 - dot(v_i_r - v_i, v_i - v_i_l)/((norm(v_i_r - v_i)*norm(v_i - v_i_l)) + eps))/2 - 1);
                            %angle between the two vectors
                            tE = [tE; left, curv{left, conf(num) - 1}(1, 2) + alpha*Edist + beta*Esmo + kappa*Eext, right];
                        end
                    end
                    minpos = find(min(tE(:, 2)) == tE(:, 2));
                    curv{k, conf(num)} = tE(minpos(1), :);
                end
            else
                %handling the last point
                curv{1, snakePt} = [1, curv{1, snakePt - 1}(1, 2) + kappa*normE(1, snakePt), 1];
            end
        end
    end

    snakeIx = []; snakeIy = [];

    minEn = 100;
    for k = 1:size_a(numP);
        if curv{k, numP}(1, 2) < minEn
            minEn = curv{k, numP}(1, 2); minpos = k;
        end
    end
    
    snakeIx = [posx(minpos, numP), snakeIx];
    snakeIy = [posy(minpos, numP), snakeIy];

    %modifying position
    for snakePt = numP:-1:2
        snakeIx = [posx(curv{minpos, snakePt}(1, 1), snakePt - 1), snakeIx];
        snakeIy = [posy(curv{minpos, snakePt}(1, 1), snakePt - 1), snakeIy];
        minpos = curv{minpos, snakePt}(1, 1);
    end

    [krigA, krigB] = krige(snakeIx, medfilt1(snakeIy), 1, size(I, 2));

    % using a least sqaures spline
    spl2 = spap2(knots, 4, krigA', krigB');
    x = 1:size(I, 2);
    y = spval(spl2, x);
    
    % increasing th eknots for the B-spline
    knots = knots * 2;    

    % increasing the points for the snake
    subsample = subsample/2;
    % decreasing search space assumig localized contour
    probWy = probWy/2;

    %median filtering possible outliers
    snakeIx = x(1:subsample:end);
    snakeIy = y(1:subsample:end);
end

% resample the data and update the paramters
a = linComb;
[xx, yy] = meshgrid(1:size(I, 2), 1:size(I, 1));
t = yy < repmat(y, size(yy, 1), 1);
a(t) = NaN;
%estimating the mixture of distributions
[f, xi] = ksdensity(a(:));
%using spline to interpolate the pdf
sp = spapi(5, xi, f);

xx = x; yy = y;