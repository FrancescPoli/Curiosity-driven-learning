function [outcome, cp, distMean]=outcomegen(numOutcomes,sigma,mean,Haz,safe,valrange)

% pick some 
% how long should the block of trials be?
% standard deviation of the generative dist...
% probability of a change-point on any given trial
% except that we set hazard rate equal to zero for "safe" trials after a change-point

% drift: if zero, no drift, if >0 it indicates the strength of the drift


% generate outcomes
maxval=valrange(end);
%mean=round(rand(1).*maxval);
outcome=NaN(numOutcomes, 1); % this will be an array of outcomes
distMean=NaN(numOutcomes, 1);% this will be an array of distribution mean
cp=zeros(numOutcomes, 1);     % this will be an array of binary change-point variable
s=safe;


for i = 1:numOutcomes
    if rand<Haz && s==0
        mean=round(rand(1).*maxval);
        cp(i)=1;
        s=safe;
    else
        s=max([s-1, 0]);
    end
        while ~isfinite(outcome(i))||outcome(i)>maxval||outcome(i)<1
            outcome(i)=round(normrnd(mean, sigma));
        end
        distMean(i)=mean;
end

end



