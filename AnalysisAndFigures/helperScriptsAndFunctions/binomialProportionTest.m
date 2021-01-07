function [pValues,zStats,hValues] =...
    binomialProportionTest(eventCount,totalCount)

% This function will calculate the z statistic from the binomial proportion 
% test for the success rates between each reward/punishment condition. It 
% uses the two-tailed binomial proportion test. The z stat needs to be 
% greater than 1.96 or less than -1.96 to have p < 0.05. The function will
% also calculate the p values for each z statistic using the normcdf
% function. 
% A good overview of the binomial proprtion test is here (as of 7/24/20):
%
% https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/binotest.htm
%
% The goal of the test is to compare two "proportions" (here, trial status
% rates) of a binary event (the status happened or it did not) occuring
% between two groups (reward sizes). The null hypohtesis is that the two
% proportions are equal. I think it implicitly assumes that the proportion
% can be seen as the Bernoulli parameter for the given group.
%
% At a glance, I got similar results for this vs. a two-tailed
% bootstrapping test (this is a bit higher power), so it checks out.
%
% correctTrials = 1 x numRewards vector of the number of correct trials to
% each reward
%
% numberTimesSeen = 1 x numRewards vector of the total number of trials to
% each reward
%
% Written by Nick Pavlovsky, validated and edited by Adam Smoulder, 7/24/20

% Only take row vectors!
eventCountSize = size(eventCount);
totalCountSize = size(totalCount);
if numel(eventCountSize) > max(eventCountSize) || numel(totalCountSize) > max(totalCountSize)
    error('Inputs are too big - should be just a row vector!')
end

% ...but if you pass in a column vector, we fix it...
if eventCountSize(1) > eventCountSize(2)
    eventCount = eventCount'; 
    eventCountSize = flip(eventCountSize);
end
if totalCountSize(1) > totalCountSize(2)
    totalCount = totalCount'; 
    totalCountSize = flip(totalCountSize); % just for bookkeeping
end

% Make a matrix to store the z statistic in. 
zStats = NaN(eventCountSize(2),eventCountSize(2));
pValues = NaN(size(eventCount,2),eventCountSize(2));
hValues = NaN(eventCountSize(2),eventCountSize(2));

% Run the binomial proportion test on the counts for each condition.
for i = 1:eventCountSize(2)
    n1 = totalCount(1,i);
    p1 = eventCount(1,i)/n1;
    for j = 1:size(eventCount,2)
        n2 = totalCount(1,j);
        p2 = eventCount(1,j)/n2;
        phat = (n1*p1 + n2*p2)/(n1 + n2);
        zStats(i,j) = (p1 - p2)/sqrt(phat*(1-phat)*(1/n1 + 1/n2));
        
        % Calculate the p-value for this z statistic. Need to multiply by 2
        % since this is a two-tailed test (H0: p1 = p2).
        pValues(i,j) = 2*(1 - normcdf(abs(zStats(i,j))));
        
        if abs(zStats(i,j)) >= 1.96
            hValues(i,j) = 1;
        else
            hValues(i,j) = 0;
        end
    end
end

end