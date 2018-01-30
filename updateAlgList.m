function updateAlgList
global algs NumOfMode;
global tdalabStatus;
TDModelstr=defstr('TDModel');
LA=length(algs);
activeAlgs=ones(1,LA);

hShowAllAlgs=findobj('style','checkbox','tag','cbShowAllAlgs');
if ~get(hShowAllAlgs,'value')
    if ~isempty(tdalabStatus.model)
        n=find(strcmpi(TDModelstr,tdalabStatus.model)==1);
        activeAlgs=activeAlgs&([algs(:).model]==n);
    end

    if ~isempty(tdalabStatus.inputType)
        activeAlgs=activeAlgs&([algs(:).numOfMode]>=max(NumOfMode));
    end

    if tdalabStatus.nonnegativity
        activeAlgs=activeAlgs&([algs(:).nonnegativity]~=AlgType.No);
    else
        activeAlgs=activeAlgs&([algs(:).nonnegativity]~=AlgType.Yes);
    end
else
    activeAlgs=activeAlgs==1;
end

algstr={algs(activeAlgs).details};
flag=1:LA;
flag(~activeAlgs)=[];
h=findobj('tag','pmAlgList');
if ~isempty(flag)
    algIndex=flag(1);
    set(h,'string',algstr,'value',1);
    h=findobj('tag','pbOpts');
    set(h,'enable','on');
else
    algstr={'No valid algorithms available'};
    algIndex=[];
    set(h,'string',algstr,'value',1);    
    h=findobj('tag','pbOpts');
    set(h,'enable','off');
end
tdalabStatus.algIndex=algIndex;
tdalabStatus.validAlgs=flag;
end