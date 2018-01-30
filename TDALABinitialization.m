function TDALABinitialization
global numOfMode numOfComp PMFalgIDs PMFalgParaList;
global Y Ycap SIRs MCTimes;
global algOptions algs;
global tdalabStatus;

algs=[];
algsInitialize;
algOptions=cell(1,numel(algs));

PMFalgIDs=[];PMFalgParaList=[];

Y=[];
Ycap=[];SIRs=[];
numOfComp=inf;numOfMode=0;

Ynoise=[];SNR=inf;
global TEMPFILE;
save(TEMPFILE,'Ynoise','SNR');

MCTimes=0;

tdalabStatus=struct('inputType',[],... % original input type
    'fullTensor',false,...  % real type for decomposition
    'model',[],...  % CP or Tucker
    'nonnegativity',false,...
    'algIndex',[],...
    'validAlgs',1:numel(algs),...
    'noiseType',[],...
    'noiseSNR',inf,...
    'noiseSparsity',0,...
    'decomposed',false,...
    'advEvaluation',false);

%% UI initialize
set(findobj('tag','chkmoderank'),'value',false);
