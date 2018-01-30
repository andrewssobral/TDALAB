function hmainOpts=setOpts(algIndex)
global tdalabStatus  paraList algOptions paraTypeList algs;

global NumOfComp;


if nargin==0
    tdalabStatus.algIndex=tdalabStatus.algIndex(1);
    algIndex=tdalabStatus.algIndex;
end
if isempty(algIndex)
    errordlg('No algorithm selected.','Error','modal');
    return;
end


for idx=algIndex
    if ~isstruct(algOptions{idx})
        algOptions{idx}=struct();
    end
    
    if isfield(paraList{idx},'NumOfComp')&&~isempty(NumOfComp)&&sum(NumOfComp)>0
        algOptions{idx}.NumOfComp=NumOfComp;
    end
end

[opts hmainOpts]=guiSetOpts(algs(algIndex),paraTypeList(algIndex),paraList(algIndex),algOptions(algIndex));

if ~isempty(opts)
    pos=1;
    for idx=algIndex
        algOptions{idx}=opts{pos};
        pos=pos+1;
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function setPMFalgOpts(hObject,eventdata,source)
%             flag=ids>1;
%             actModes=1:NumOfMode;
%             actModes=actModes(flag);
%             if sum(flag)>0
%                 PMFActParas=PMFsetOpts((ids(flag))-1,actModes);
%             end
%         end
%         function exitBSSAlgID(hObject,eventdata,source)
%             ids=ids(:);
%             set(heditBSSAlgID,'string',num2str(ids'));
%             
%            %% set default options for PMF here if it is empty
%             if isempty(PMFActParas)
%                 flag=ids>1;
%                 PMFActParas=PMFParaListParsing((ids(flag))-1);
%             end
%             
%             delete(hfigAlgID);
%         end   