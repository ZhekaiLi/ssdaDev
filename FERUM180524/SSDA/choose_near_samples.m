function [sampleX, sampleY] = choose_near_samples(baseSampleX, baseSampleY, newSampleX, newSampleY)

    n = size(baseSampleX, 2);
    m = size(newSampleX, 2);

    meanX = mean(newSampleX, 2);

    Dis = zeros(1, n);

    for i = 1:n
        dis = norm(baseSampleX(:,1) - meanX);
        Dis(1, i) = dis;
    end

    q = quantile(Dis, (n-m)/n);
    index =  find(Dis <= q);
    newX = baseSampleX(:, index);
    newY = baseSampleY(:, index);

    sampleX = [newX, newSampleX];
    sampleY = [newY, newSampleY];
end