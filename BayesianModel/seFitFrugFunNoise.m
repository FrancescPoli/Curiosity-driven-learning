function [estimates, modPred, sse, estimateserror, totSig, bic, R, pCha, Alpha] = seFitFrugFunNoise (xdata, ydata, nn, haz, whichParams, tR, zeroNode, newBlock)
% whichParams is a logical array, length=4, specifying which parameters
% you want to fit.  1 = hazard, 2= noise, 3=drift, 4=likelihood weight.
%
% xdata = outcomes
% ydata = subject estimates
% conf = subject confidence
% windows (if measured... leave empty if not used)
% nn = standard deviation of distribution
% haz = hazard rate
% tR = trueRun (see frugFun5... basically determines exactly how variance is updated by the model)
% zeroNode = 1 includes uniform prior distribution in all belief estimates (0 for better fit of subject data).
% newBlock = binary vector, true for each trial that was the first in a new block, false everywhere else


wp=logical(whichParams);
sp=[haz, nn(1),   0,   1];
lb=[ 0 , .001, 0,   0];
ub=[ 1 , 100, 100,  1];

if nargin < 8 | isempty(zeroNode)
    zeroNode=0;
end

if nargin < 7 | isempty(tR)
    tR=0;
end

% so stupid.  don't have time to do a better job.
if size(ydata, 2)==1
    ydata=ydata'
end


% for some reason fmincon gives up when fitting the drift if it starts at
% zero.  lets just start it at something bigger.
if wp(3)
    sp(3)=20;
end


%% Call fmincon to evaluate squred error in fitting a curve with a parameters initially set to starting point.

start_point = sp(wp)
model = @frugFun;                % this is the function handle to the function that takes the parameters and outputs the thing we want to minimize

 oldOpts = optimset('fmincon');
 options=optimset(oldOpts, 'maxFunEvals', 1000000000000, 'MaxIter', 10e10);
[estimates, sse, ef, o, l, g, h] = fmincon(model, start_point, [], [], [], [], lb(wp), ub(wp), [], options);
% this calculates the error from the hessian
estimateserror = sqrt(diag(-((-h)^(-1))));

%% now do get info for *BEST* fit model
begs=find(newBlock);
ends=[begs(2:end)-1 length(newBlock)];
goodPar(wp) = estimates;
goodPar(~wp)= sp(~wp);
goodPar(1)  = max([goodPar(1), 0]);
for i = 1:length(begs)
    if wp(2)
        noise=goodPar(2);
    else
        noise=nn(i);
    end
    sel=begs(i):ends(i);
    [B, C, D, E, ~, F] = frugFunNoise(xdata(sel), goodPar(1), noise, goodPar(3), goodPar(4), 1, length(newBlock),tR, ydata(begs(i)));
    modPred(begs(i):ends(i)) = B(1:ends(i)-begs(i)+1);                  %function mapping inputs to estimates
    totSig(begs(i):ends(i)) = C(1:ends(i)-begs(i)+1);
    R(begs(i):ends(i))=D(1:ends(i)-begs(i)+1);
    pCha(begs(i):ends(i))=E(1:ends(i)-begs(i)+1);
    Alpha(begs(i):ends(i))=F(1:ends(i)-begs(i)+1);
end

ErrorVector = modPred - ydata;                     %error at each point
ErrorVector(~isfinite(ErrorVector))=300;
sse = nansum(ErrorVector.^2);
% bayesian information criterion... applies stiff punishment for more fit
% terms:
bic=length(ydata).*log(sum(ErrorVector.^2)./length(ydata))+sum(whichParams).*log(length(ydata));




% frugFun accepts some inputs (params) and an array that tells it which
% inputs they are (wp)
    function [totDist] = frugFun(params)                %input comes from fmincon (initially start_point)
        
        begs=find(newBlock);
        ends=[begs(2:end)-1 length(newBlock)];
        
        
        goodPar(wp) = params;
        goodPar(~wp)= sp(~wp);
        goodPar(1)  = max([goodPar(1), 0]);
         
        for i = 1:length(begs)
            if wp(2)
                noise=params(2);
            else
                noise=nn(i);
            end 
            sel=begs(i):ends(i);
            %keyboard
            [B, totSig,~,~,~,~] = frugFunNoise(xdata(sel), goodPar(1), noise, goodPar(3), goodPar(4), 1, length(newBlock),tR, ydata(begs(i)));
            
            modPred(begs(i):ends(i)) = B(1:ends(i)-begs(i)+1);                  %function mapping inputs to estimates
        end
         
        mu1 =modPred;
        mu2 =ydata;
        
        
        
%         
%         sqErr=(mu1-mu2).^2;
%         sqErr(~isfinite(sqErr))=300.^2;
%         totDist = nansum(sqErr)        
%     
        absErr=abs(mu1-mu2);
        absErr(~isfinite(absErr))=300;
        totDist=nansum(absErr);
    
    end

end




