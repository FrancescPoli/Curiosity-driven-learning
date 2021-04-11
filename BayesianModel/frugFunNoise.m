function [B, totSig, R, pCha, vari, Alpha] = frugFunNoise(data, hazExp, noise, drift, likeWeight, learnNoise, C, trueRun, initGuess)  
%input comes from fmincon (initially start_point)
% data is the values to estimate, hazExp is the expected hazard rate, noise is the expected standard deviation of the
% distribution, drift is the standard deviation of a mean-centered drift
% distribution, likeWeight assigns a weight to the likelihood relative to
% the prior


if nargin <9 | isempty(initGuess)
    initGuess=22;
end

if hazExp <0
    hazExp = 0;
end

if hazExp >1
    hazExp =1;
end

Alpha=[];
sigmaE = noise;
vari    = sigmaE.^2;
B(1)   = initGuess;
R(1)   = 1;
nCount = 0;
%expRun=.1;

% 
% 
% % determine A and B parameters of beta distribution based on expected value
% % and sharpness
% HazA=hazExp.*hazSharp;
% HazB=hazSharp-HazA;


if nargin < 5 | isempty(likeWeight)
    likeWeight=1;
end

if nargin < 4 | isempty(drift)
    drift=0;
end



%% this is stupid code to get the "true prior"   This stuff should really go inside 
% a=zeros(300,1);                                 % where means were picked...
% a(71:230)=1    ;                                % where means were picked...
% b=normpdf(1:301, 151, noise);                   % uncertainty due to noise
% c=conv(a, b)   ;                                % likelihood of getting a number given that a changepoint occurred
% d=c(ceil(length(b)./2):length(a)+floor(length(b)./2));
% d=d./sum(d);
% 
%% this is a flat prior...
d=ones(43,1)./43;



clear chaRat
%% loop through data making sequential predictions

for i = 1:length(data)
    
    if learnNoise==1
    % reset sigmaE based on errors    
    sigmaE=sqrt(vari(i));
    end


    % part 1 get the expected distribution
    sigmaU(i)=sqrt(((sigmaE.^2)./R(i))+drift.^2);         % changed things to squared and added in drift uncertainty
    R(i)=sigmaE.^2 ./sigmaU(i).^2;                       % recompute R including drift uncertainty
    totSig(i)=sqrt(sigmaE.^2+sigmaU(i).^2);             % same deal
     

    % part 2 calculate probability of change
    pI=normpdf(data(i),B(i),totSig(i));
    % normalize to correct for probability outside of range.  fixed
    % an error in this code on 9-25-09...
    normalize=((normcdf(43, B(i), totSig(i)))- normcdf(0, B(i), totSig(i)));
    pI=pI./normalize;


    if data(i) < 1 || data(i) > 43                      % compute likelihood of getting data given that changepoint occurred
        changLike=d(1);
    else
        changLike=d(data(i));
    end

    %changeRatio=exp(likeWeight.*log(changLike./pI)+ (log(H1)-log(Hmin1)));
    
    likeRatio  =log(changLike)-log(pI);  
    prioRatio  =log(hazExp)-log(1-hazExp);
    changeRatio= exp(likeWeight.*likeRatio + prioRatio);
    chaRat(isfinite(changeRatio))=changeRatio(isfinite(changeRatio));
    chaRat(~isfinite(changeRatio))=exp(100);
    allPCha=chaRat'./(1+chaRat');
    pCha(i)=nanmean(allPCha);
    pNoCha=1-pCha(i);

    % part 3 update belief about mean

    yInt  = 1./(R(i)+1);    % now R can be really small if there is a big drift...
    slope = (1-yInt);
    Alph  = yInt+pCha(i).*slope;
    Alpha(i)=Alph;
    Delta = data(i)-B(i);
    B(i+1)= B(i)+Alph.*Delta;

    % part 4 update run length   
    if trueRun==1
        R(i+1)=(R(i)+1).*pNoCha+pCha(i);% computed by taking expected value over the run length distribution
    else
        ss=pCha(i).*((sigmaE.^2)./1)+pNoCha.*((sigmaE.^2)./(R(i)+1))+pCha(i).*pNoCha.*((B(i)+yInt.*Delta)-data(i)).^2;
        R(i+1)=(sigmaE.^2)./ss;         % or computed by matching moments of predictive distribution
    end
   
    % part 5 update estimate of hazard rate

%    if learnHaz==1
        % to use this code use frugFun8 script.  Though it may not be very
        % up to date... im just sayin...
%     I       =HazA+HazB+1;
%     mu1     =(HazA+1)./(I);
%     mu2     =(HazA)./(I);
%  
%     var1    =mu1.*(1-mu1)./(I+1);
%     var2    =mu2.*(1-mu2)./(I+1);
%     var3    =(mu1-mu2).^2;
%     
%     muE     =mu1.*pCha(i)+mu2.*pNoCha;
%     varE    =pCha(i).*var1+pNoCha.*var2+pCha(i).*pNoCha.*var3;
%     
% 
%     HazA=muE.*(muE.*(1-muE)-varE)./varE;
%     HazB=(1-muE).*(muE.*(1-muE)-varE)./varE;
%   
%     % add maximum certainty to hazard rate
%     HazA     =max([HazA-(HazA./C) 0]);
%     HazB     =max([HazB-(HazB./C) 0]);
%     expHaz(i)=HazA./(HazB+HazA);
%     sharp(i) =HazA+HazB;
%     
%     else
%         expHaz(i)=hazExp;
%         sharp(i) =hazSharp;
%     end
    
    
    % part 6: update estimate of noise  (THIS SHOULD BE UNBIASED!!!)    
    
    if i >1
    exVar   = (R(i)./(R(i)+1)).*(Delta.^2);
    nCount  = nCount+(1-pCha(i));   
    vari(i+1)= (vari(i).*(nCount-pNoCha)./nCount)       +       exVar.*pNoCha./nCount;  % this is an updated assessment of the variance
    
    %nCount=max([nCount-nCount./C 0]);
    
    %nCount=max([0 min([nCount 1./C-1])]);
    
    else 
        vari(i+1)=vari(i);
    end
  
    
   
end


end