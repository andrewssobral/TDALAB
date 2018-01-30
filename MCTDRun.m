function MCTDRun
global algOptions algs tdalabStatus paraList hmainTDALAB ;
global Y Ycap NumOfComp NumOfMode SIRs MCfit MCSIRs MCTimes MCelapsedTime;
TDModelstr=defstr('TDModel');
commandwindow;
if isempty(tdalabStatus.algIndex)
    errordlg('Please select an algorithm first.','Setting error','modal');
    return;
end

fprintf('[TDALAB] Initializing for tensor decomposition ... \n');
%% initialize
SIRs=[];Ycap=[];Ynoise=[];

fprintf('[TDALAB] Checking the parameters ...\n');
%% options
for algIndex=tdalabStatus.algIndex;
    if isempty(algOptions{algIndex})
        algOptions{algIndex}=struct;
    end

    if isfield(paraList{algIndex},'NumOfComp')&&(~isfield(algOptions{algIndex},'NumOfComp'))
        algOptions{algIndex}.NumOfComp=NumOfComp;
    end
    if isfield(algOptions{algIndex},'NumOfComp')&&strcmpi(tdalabStatus.model,'CP')
        algOptions{algIndex}.NumOfComp=algOptions{algIndex}.NumOfComp(1);
    end

    fnames=fieldnames(algOptions{algIndex});
    for n=1:numel(fnames)
        if isempty(algOptions{algIndex}.(fnames{n}))
            algOptions{algIndex}=rmfield(algOptions{algIndex},fnames{n});
        end
    end
end

%% load observations Ynoise [tdalabtemp.mat]
global TEMPFILE;
fprintf('[TDALAB] Preparing the data for decomposition ... \n');
if exist(TEMPFILE,'file')~=0
    load(TEMPFILE,'Ynoise','SNR');
else
    Ynoise=[];SNR=rand;
end
if isempty(Ynoise)||(SNR~=tdalabStatus.noiseSNR)
    if tdalabStatus.fullTensor
        if isinf(tdalabStatus.noiseSNR)
            fprintf('[TDALAB] Generating the full tensor. This may cost a few minutes depending on the problem size. Please wait ...\n');
            Ynoise=tensor(Y);
            SNR=tdalabStatus.noiseSNR;
            save(TEMPFILE,'Ynoise','-v7.3','SNR');
        else
            errordlg('Noise data not found. Please check your noise settings.','Invalid input','modal');
            return;
        end
    else
        Ynoise=tensor(Y);
    end
end
%% disp basic information
if strcmp(tdalabStatus.model,'CP')
    fprintf('===============   Tensor decomposition [CP model] =============== \n');
else
    fprintf('===============  Tensor decomposition [Tucker model] =============== \n');
end

fprintf(horzcat('* Input type           : ',tdalabStatus.inputType,'\n'));
str=horzcat('* Dimension of tensorb : %d',repmat(' x %d',1,NumOfMode-1),'\n');
fprintf(str,size(Y));
str=horzcat('* Number of Components : %d',repmat(' x %d',1,NumOfMode-1),'\n');

fprintf('Compared algorithms:\n');
i=1;
for algIndex=tdalabStatus.algIndex    
    dname=algs(algIndex).details;
    dname=regexp(dname,'</.*>','split');
    dname=regexp(dname{1},'<.*>','split');
    dname(cellfun(@isempty,dname))=[];
    fprintf('%d. ------------------------------------------------------\n',i);
    fprintf(horzcat('* Algorithm : ',dname{1},'\n'));
    i=i+1;
    fprintf('* Options:\n');
    disp(algOptions{algIndex});
end

%% begin to run
tdalab('hide');
NAlgs=numel(tdalabStatus.algIndex);
MCSIRs=NaN(NAlgs,NumOfMode,MCTimes);
MCfit=zeros(NAlgs,MCTimes);
MCelapsedTime=zeros(NAlgs,MCTimes);
allR=NAlgs*MCTimes;
i=1;
done=1;
for algIndex=tdalabStatus.algIndex
    currentAlg=str2func(algs(algIndex).name);
    for mc=1:MCTimes        
        fprintf('[TDALAB] %d%% completed. Alg[%s]: %d/%d. Repeat: %d/%d.\n',round((done-1)*100/allR),algs(algIndex).name,i,NAlgs,mc,MCTimes);
        ts=tic;
        Ycap=currentAlg(Ynoise,algOptions{algIndex});
        MCelapsedTime(i,mc)=toc(ts);
        if ~isempty(Ycap)&&(~strcmpi(tdalabStatus.inputType,'tensor'))
            for n=1:NumOfMode
                if size(Y.U{n},2)==size(Ycap.U{n},2)
                    MCSIRs(i,n,mc)=mean(CalcSIR(Y.U{n},Ycap.U{n}));
                else
                    MCSIRs(i,n,mc)=NaN;
                end
            end
        end
        MCfit(i,mc)=fitness(Y,Ycap);
        done=done+1;
    end
    i=i+1;
end

set(hmainTDALAB,'HandleVisibility','callback','visible','on');
tdalabStatus.decomposed=true;
commandwindow;

MCVisualize;

fprintf('\n========================== END ==========================\n\n');
fprintf('\n[TDALAB] Type or click <a href="matlab: tdalab">tdalab</a> to return.\n\n');

if tdalabStatus.decomposed
    set(allchild(findobj('tag','pnlOutputAnalysis')),'enable','on');
else
    set(allchild(findobj('tag','pnlOutputAnalysis')),'enable','off');
end

    



end
