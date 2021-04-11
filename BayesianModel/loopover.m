%% Script that loops every sequence of every participant 
% to fit the model (frugFunNoise.m) to the data (matdata)

%% Import data from text file.
% Initialize variables.
filename = 'U:\Desktop\S_AdultCuriosity\BayesianModel\matdata.csv';
delimiter = ',';
startRow = 2;

%% Format for each line of text:
%   column2: double (%f)
%	column3: double (%f)
%   column4: double (%f)
%	column5: double (%f)
%   column6: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Create output variable
matdata = [dataArray{1:end-1}];
%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% Now run sefitfrugfun for every block of every participant
% Initialize variables
allsig=[];
allalpha=[];
allB=[];
alldat=[];
allsubj=[];
allseq=[];
n=0;
for i = 1:nnz(matdata(:,3)==1)
    % three sequences per subject
    n = n+1;
    if n==4
        n=1; %after sequence 3 go back to 1
    end
    idx=find(matdata(:,3)==1,i); % get data for this sequence
    idx=idx(end);
    % loop to adjust indexing
    if i < nnz(matdata(:,3)==1)
        idxend=find(matdata(:,3)==1,i+1)-1;
        idxend=idxend(end);
    elseif i == nnz(matdata(:,3)==1)
        idxend=length(matdata);
    end
    xdata = matdata(idx:idxend,2); % specify x and y data (real positions and participants' guesses)
    ydata = matdata(idx:idxend,1);
    newBlock = matdata(idx:idxend,3);
    %initialize some variables for the model
    nn = matdata(idx,4);
    haz = matdata(idx,5);
    whichParams=[1 1 0 0]; %which parameters to fit
    tR=0;
    zeroNode=0;
    % fit the model
    [estimates, modPred, sse, estimateserror, totSig, bic, R, pCha, alpha] = seFitFrugFunNoise(xdata, ydata, nn, haz, whichParams, tR, zeroNode, newBlock);
    % update all variables
    allB=[allB,modPred];
    allalpha=[allalpha,alpha];
    allsig=[allsig,totSig];
    alldat=[alldat,ydata'];
    allsubj=[allsubj,[ones(length(ydata),1)*i]'];
    allseq=[allseq,((n-1)*35+1):((n-1)*35+length(xdata))'];
end

% export model parameters
expmat=[allB;allsig;allalpha;alldat;allsubj;allseq]';
dlmwrite('expmat.csv',expmat)

% average alpha values
alphaavg=splitapply(@mean,expmat(:,3),findgroups(expmat(:,6)));
dlmwrite('alphaavg.csv',alphaavg)

