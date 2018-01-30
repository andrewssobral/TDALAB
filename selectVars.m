function selectVars
global Y NumOfMode NumOfComp loadedVars tdalabStatus loadedFileName;
TDModelstr=defstr('tdmodel');
h=findobj('tag','ppSelVar');
index=get(h,'value');
names=loadedVars;
selvar=load(loadedFileName,names{index});
if isempty(selvar)
    error('[TDALAB] The source file has been modified. Please reload it.');
end
selvar=selvar.(names{index});
clsvar=class(selvar);
hmodel=findobj('tag','pmTDModel');

inputFormatstr=defstr('input');

if ~any(strcmpi(inputFormatstr,clsvar))&&~any(strcmpi({'double','tensor'},clsvar))
    h=findobj('tag','txtNoTensorInfor');
    set(h,'visible','on');
    set(h,'string',{selvar});
elseif strcmp(clsvar,'double')||strcmp(clsvar,'tensor')
    Y=selvar;Y=tensor(Y);
    cleartemp;
    tdalabStatus.inputType='tensor';
    tdalabStatus.model='CP';
    tdalabStatus.fullTensor=true;
    sz=size(Y);NumOfMode=length(sz);
    NumOfComp=0;
    NumOfComp=NumOfComp(ones(1,NumOfMode));
    set(hmodel,'value',find(strcmpi(TDModelstr,tdalabStatus.model),1));
    updateUI;
elseif strcmp(clsvar,'ktensor')
    Y=selvar;
    cleartemp;
    tdalabStatus.inputType='ktensor';
    tdalabStatus.model='CP';
    tdalabStatus.fullTensor=false;
    NumOfMode=length(size(Y));
    NumOfComp=length(Y.lambda);
    NumOfComp=NumOfComp(ones(1,NumOfMode));
    set(hmodel,'value',find(strcmpi(TDModelstr,tdalabStatus.model),1));
    updateUI;
elseif strcmp(clsvar,'ttensor')
    Y=selvar;
    Y=ttensor(tensor(Y.core),Y.U);
    cleartemp;
    tdalabStatus.inputType='ttensor';
    tdalabStatus.model='Tucker';
    tdalabStatus.fullTensor=false;
    NumOfComp=size(Y.core);
    NumOfMode=length(size(Y));
    set(hmodel,'value',find(strcmpi(TDModelstr,tdalabStatus.model),1));
    updateUI;
else
    errordlg(['Unsuported data type: [',class(Y),'].'],'Type error','modal');
    return;
end
end