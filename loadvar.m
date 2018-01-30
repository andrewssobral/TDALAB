function v=loadvar(type,title,defaultDir,defaultVarName)
if nargin<3
    defaultDir=[];
    defaultVarName=[];
elseif nargin<4
    defaultVarName=[];
end

bgColor=[0.8314    0.8157    0.7843];
hdtype=[];
hdsize=[];
currVar=1;
vars=[];
v=[];

[FileName FilePath]=uigetfile(type,title,defaultDir);
if ~isequal(FileName,0)
    FileName=horzcat(FilePath,FileName);
    vars=whos('-file',FileName);
    
    if ~isempty(defaultVarName)
        currVar=find(strcmpi({vars(:).name},defaultVarName)==1);
    end
    
    if numel(vars)>1
        figvar=figure('Units','Characters','Resize','on','toolbar','none','menu','none','color',bgColor,...
            'OuterPosition',[0 0 80 15],...
            'tag','figTensorClss','Name','Select a variable to load ...','numbertitle','off','visible','off');
        hspace=1/5;hctrl=hspace*.75;
        % label
        uicontrol('parent',figvar,'units','normalized','position',[0.01 3.5*hspace 0.4 hctrl],'style','text',...
            'string','Select a variable:','horizontalalign','right','backgroundcolor',bgColor,...
            'fontunits','normalized','fontsize',.6);
        uicontrol('parent',figvar,'units','normalized','position',[0.01 2.5*hspace 0.4 hctrl],'style','text',...
            'string','Data type:','horizontalalign','right','backgroundcolor',bgColor,...
            'fontunits','normalized','fontsize',.6);
        uicontrol('parent',figvar,'units','normalized','position',[0.01 1.5*hspace 0.4 hctrl],'style','text',...
            'string','Date size:','horizontalalign','right','backgroundcolor',bgColor,...
            'fontunits','normalized','fontsize',.6);
        
        % infor
        uicontrol('parent',figvar,'units','normalized','position',[0.42 3.5*hspace 0.4 hctrl],'style','popupmenu',...
            'string',{vars(:).name},'callback',@selVar,...
            'fontunits','normalized','fontsize',.6,'value',currVar);
        hdtype=uicontrol('parent',figvar,'units','normalized','position',[0.42 2.5*hspace 0.4 hctrl],'style','text',...
            'string',vars(currVar).class,'horizontalalign','left','backgroundcolor',bgColor,...
            'fontunits','normalized','fontsize',.6);
        hdsize=uicontrol('parent',figvar,'units','normalized','position',[0.42 1.5*hspace 0.4 hctrl],'style','text',...
            'string', strsize(vars(currVar).size),'horizontalalign','left','backgroundcolor',bgColor,...
            'fontunits','normalized','fontsize',.6);
        
        uicontrol('parent',figvar,'units','normalized','position',[0.2 0.2*hspace 0.25 hctrl],'style','pushbutton',...
            'string','OK',...
            'fontunits','normalized','fontsize',.6,'callback',@loadcurrVar);        
        uicontrol('parent',figvar,'units','normalized','position',[0.55 0.2*hspace 0.25 hctrl],'style','pushbutton',...
            'string','Cancel',...
            'fontunits','normalized','fontsize',.6,'callback',@cancel);
        
        movegui(figvar,'center');
        set(figvar,'visible','on');
        waitfor(figvar);
        
    else
        v=load(FileName,vars(1).name);
        v=v.(vars(1).name);
    end
else
    v=[];
end

    function str=strsize(sz)
        str=regexprep(num2str(sz),'\s+',' x ');
    end
    function selVar(hObject,eventdata)
        currVar=get(hObject,'value');
        set(hdtype,'string',vars(currVar).class);
        set(hdsize,'string',strsize(vars(currVar).size));
    end
    function loadcurrVar(hObject,eventdata,handles)
        v=load(FileName,vars(currVar).name);        
        v=v.(vars(currVar).name);
        close(figvar);
    end
    function cancel(hObject,eventdata,handles)
        v=[];
        close(figvar);
    end

end