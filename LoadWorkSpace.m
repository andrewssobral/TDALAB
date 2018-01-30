function LoadWorkSpace
global TDALABHOME Y tdalabStatus;
% if isemtpy(Y)
%     errordlg('Did not load any data yet.');
% end
global NumOfComp NumOfMode;
global algNames algParas Ynoise noiseSNR;
global MCTimes MCelapsedTime MCSIRs MCfit;
global elapsedTime SIRs fit Ycap;
global loadedVars loadedFileName;
global algs algOptions;

[WSFileName WSFilePath]=uigetfile('*.mat','Load a workspace ...',horzcat(TDALABHOME,filesep,'userdata'));
if isequal(WSFileName,0)
    return;
end
MCFileName=horzcat(WSFilePath,WSFileName);
ws=load(MCFileName);

loadedFileName=ws.Exp.dataFile;
fvars=whos('-file',loadedFileName);
loadedVars={fvars(:).name};

FileName=regexp(ws.Exp.dataFile,filesep,'split');
FileName=FileName{end};
h=allchild(0);
oo = findall(h, 'Tag', 'TDALAB' );
set(oo,'name',[FileName,' :: ',defstr('TDALAB')]);

TDALABinitialization;
tdalabStatus=struct();
algParas=struct();

fs=fieldnames(ws.Exp);
for idx=1:numel(fs)
    eval([fs{idx} '=ws.Exp.(fs{idx});']);
end

%% load source tensor data
data=load(ws.Exp.dataFile);
if isempty(data)
    errstr=horzcat('The file ',ws.Exp.dataFile,' does not exist.');
    errordlg(errstr,'Open error');
    return;
end
varnames=fieldnames(data);
idx=find(strcmpi(varnames,ws.Exp.Y));
h=findobj('tag','ppSelVar');
set(h,'string',varnames,'enable','on','value',idx);
Y=data.(varnames{idx});

NumOfMode=numel(size(Y));
if ~isempty(Ycap)
    switch lower(tdalabStatus.model)
        case {'cp','tucker'}
            NumOfComp=reshape(cellfun(@(x) size(x,2),Ycap.U(:)),1,[]);
        case {'bcd'}
    end
            
else
    if find(strcmpi({'ttensor','ktensor'},class(Y)))
        NumOfComp=reshape(cellfun(@(x) size(x,2),Y.U(:)),1,[]);
    else
        NumOfComp=[];
    end
end

updateUI;


%% update algs information
if ~tdalabStatus.advEvaluation
    h=findobj('tag','pmAlgList');
    idx=find(strcmpi({algs(tdalabStatus.validAlgs).name},ws.Exp.algNames));
    if idx>0
        set(h,'value',idx);
        idx= strcmpi({algs(:).name},ws.Exp.algNames);
        algOptions{idx==1}=ws.Exp.algParas;
        tdalabStatus.algIndex=find(idx==1);
    else
        errordlg(horzcat('The algorithm [', ws.Exp.algNames, '] has been removed.'));
    end
else  % MC RUN
    %% update tdalabStatus.algIndex
    pos=1;
    for idx=1:numel(ws.Exp.algNames)
        cidx=find(strcmpi({algs(:).name},ws.Exp.algNames{idx}));
        if cidx>0
            tdalabStatus.algIndex(pos)=cidx;
            algOptions{cidx}=ws.Exp.algParas.(ws.Exp.algNames{idx});
            pos=pos+1;
        end
    end
end %% update algs

%%command visualization
commandwindow;
tdalab('hide');
fprintf('Algorithm: %s [%s].\n',algs(tdalabStatus.algIndex).details,algs(tdalabStatus.algIndex).name);
fprintf('Elapsed time of this run: %f sec. Fit=%3.2f%%.\n\n',elapsedTime,fit*100);

if ~all(isnan(SIRs))
     SIRs=real(SIRs);
     CPdispSIRs(SIRs);
end

fprintf('========================== Load complete ==========================\n\n');
fprintf('\n[TDALAB] Type or click <a href="matlab: tdalab">tdalab</a> to return.\n\n');
pause(1);
tdalab('show');

end