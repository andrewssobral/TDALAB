function  t=refineTucker()
%% Impose constraits on the Tucker decomposition
global Y Ycap tdalabStatus NumOfMode NumOfComp SIRs PMFalgIDs PMFActParas ;
ExtraConStr=defstr('TuckerConstr');
h=findobj('tag','gbExtraConstraints');
% PMFalgIDs=getappdata(h,'PMFalgIDs');
hsel=get(h,'SelectedObject');
consType=get(hsel,'string');
if isempty(Ycap)&&strcmpi(tdalabStatus.inputType,'ttensor')
    Ycap=Y;
end
if isempty(Ycap)||~strcmpi(class(Ycap),'ttensor')
    errordlg('A ttensor is required in Tucker refinement.','Error data type','modal');
    return;
end

consType=lower(consType);
cmodes=find(PMFalgIDs>1)';

[temp PMFalgs PMFparaTypeList]=PMFalgInit;


commandwindow;
% fprintf('## Type:%s. Modes:[%s]. Before refinement: fit=%f. \n',consType,num2str(cmodes,'%2d'),fit);

A=Ycap.U;G=Ycap.core;
t=tic;
switch consType
    case 'penalized matrix factorization'
        for i=1:numel(cmodes)
            n=cmodes(i);
            curralg=str2func(PMFalgs(PMFalgIDs(n)-1).name);
            PMFActParas{i}.NumOfComp=size(Ycap.U{n},2);
            fprintf('Mode %d is updated by using [%s]\n',n,PMFalgs(PMFalgIDs(n)-1).details);
            [A{n} w]=curralg(Ycap.U{n},PMFActParas{i});
            G=ttm(G,w,n);
        end
        
    otherwise
        fprintf(2,'[TDALAB] Unsuported method. Return without action.\n');
end
t=toc(t);
Ycap=ttensor(G,A);
[fit]=fitness(Y,Ycap);
fprintf('## Refinement complete: fit=%f\n',fit);
fprintf('\n======================= <a href="matlab: tdalab">TDALAB</a> =======================\n\n');   

if (~strcmpi(tdalabStatus.inputType,'tensor'))
    %         SIRs = cellfun(@CalcSIR,Ycap.U(:),Y.U(:),'uni',0)';
    SIRs=NaN(max(NumOfComp),NumOfMode);
    for n=1:NumOfMode
        if size(Y.U{n},2)==size(Ycap.U{n},2)
            SIRs(1:NumOfComp(n),n)=CalcSIR(Y.U{n},Ycap.U{n});
        end
    end
    if ~all(isnan(SIRs))
        SIRs=real(SIRs);
        CPdispSIRs(SIRs);
    end
end