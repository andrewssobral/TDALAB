function setExtraCons(source,eventdata)
global ModeConIndex NumOfMode tdalabStatus;
global lFontSize sFontSize defaultFontName PMFalgIDs PMFActParas;
global ScreenWidth ScreenHeight XSpread defaultCtrlHeight backGroundColor;
currentConstr=get(eventdata.NewValue,'String');

if strcmp(currentConstr,'None')
    ModeConIndex=[];
    return;
elseif isempty(tdalabStatus.inputType)
    return;
end

PMFActParas=[];

width=max(100*XSpread,(NumOfMode+2)*10*XSpread);
height=14*defaultCtrlHeight;
helpHeight=0;

hmainfig=figure('Units','Characters','position',[ScreenWidth/2-width/2 ScreenHeight/2-height/2 width height],...
    'Name','Extra constraints for Tucker decomposition','numbertitle','off','tag','TuckerConsWin',...
    'resize','on','toolbar','none','menu','none','windowstyle','modal','visible','off','color',backGroundColor);
set(hmainfig,'CloseRequestFcn',@cancelTuckerCons);
bkcolor=get(hmainfig,'color');
pos=get(hmainfig,'position');
%width=pos(3);height=pos(4);
nBorderWidth=0.01;
dnBorderWidth=0.02;

uicontrol('parent',hmainfig,'units','normalized','position',[nBorderWidth 0.8 1-dnBorderWidth 0.15],'string',horzcat('Constraint type: ',currentConstr),...
    'style','text','backgroundcolor',bkcolor,'fontsize',lFontSize+2,'fontName',defaultFontName);

%% select algorithms
uicontrol('parent',hmainfig,'units','normalized','position',[nBorderWidth 0.68 1-dnBorderWidth 0.1],'string','Select algorithms for each mode:',...
    'style','text','backgroundcolor',bkcolor,'fontsize',lFontSize,'fontName',defaultFontName,'horizontalalignment','left');
uicontrol('parent',hmainfig,'units','normalized','position',[nBorderWidth 0.4 1-dnBorderWidth 0.06],'string','*[Free] means no constraints on this mode',...
    'style','text','backgroundcolor',bkcolor,'fontsize',sFontSize+1,'fontName',defaultFontName,'horizontalalignment','left');
width=(1-dnBorderWidth)/NumOfMode;
for n=1:NumOfMode
    uicontrol('parent',hmainfig,'units','normalized','position',[nBorderWidth+(n-1)*width 0.58 width*0.8 0.08],'string',horzcat('Mode ',num2str(n)),...
        'style','text','backgroundcolor',bkcolor,'fontsize',lFontSize,'fontName',defaultFontName,'horizontalalignment','left');
    uicontrol('parent',hmainfig,'units','normalized','position',[nBorderWidth+(n-1)*width 0.48 width*0.8 0.1],'string',defstr(currentConstr),...
        'style','popupmenu','fontsize',lFontSize,'fontName',defaultFontName,'callback',@selModeAlgID,'userdata',n,'createFcn',@selModeAlgID);
end

hbpPMFopts=uicontrol('parent',hmainfig,'units','normalized','position',[nBorderWidth,0.28 0.3 0.1],'string','Algorithm Options',...
    'style','pushbutton','fontname',defaultFontName,'fontsize',lFontSize,'enable','on','callback',@setPMFOpts);
    

hcbUpdateMode=uicontrol('parent',hmainfig,'units','normalized','position',[nBorderWidth 0.15 0.98 0.1],...
    'backgroundcolor',bkcolor,'style','checkbox','tag','ckMode','string','Update current results directly',...
    'fontsize',lFontSize,'fontname',defaultFontName,'Max',true,'Min',false,'enable','off');
if tdalabStatus.decomposed||strcmpi(tdalabStatus.inputType,'ttensor')
    set(hcbUpdateMode,'enable','on');
end
uicontrol('parent',hmainfig,'units','normalized','position',[0.2 0.01 0.2 0.12],...
    'style','pushbutton','tag','ckMode','string','Confirm',...
    'fontsize',lFontSize,'fontname',defaultFontName,'callback',@setCons);
uicontrol('parent',hmainfig,'units','normalized','position',[0.6 0.01 0.2 0.12],...
    'style','pushbutton','tag','ckMode','string','Cancel',...
    'fontsize',lFontSize,'fontname',defaultFontName,'callback',@cancelTuckerCons);

movegui(hmainfig,'center');
set(hmainfig,'visible','on');

%% callbacks
    function setPMFOpts(hObject,eventdata,handles)
        flag=PMFalgIDs>1;
        if sum(flag)>0
            %% call guiSetOpts              
              [PMFparaList PMFalgs PMFparaTList]=PMFalgInit();
              algIds=PMFalgIDs(flag)-1;
              opts=guiSetOpts(PMFalgs(algIds),PMFparaTList(algIds),PMFparaList(algIds));
              pos=1;
              PMFActParas=cell(1,numel(algIds));
              for cidx=1:numel(algIds)
                  PMFActParas{cidx}=opts{pos};
                  pos=pos+1;
              end
        end
    end
    function selModeAlgID(hObject,eventdata,handles)
        PMFalgIDs(get(hObject,'userdata'))=get(hObject,'value');
    end
    function cancelTuckerCons(hObject, eventdata, handles)
        delete(hmainfig);
    end
    function setCons(hObject, eventdata, handles)
        if isempty(PMFActParas) %% if the parameters are not set yet
            flag=PMFalgIDs>1;
            indices=PMFalgIDs(flag)-1;
            [PMFparaList PMFalgs PMFparaTypeList]=PMFalgInit();
            PMFparaList=ParaListParsing(PMFparaList,PMFparaTypeList,unique(indices));
            pos=1;
            for cidx=indices
                PMFActParas{pos}=PMFparaList{cidx};
                pos=pos+1;
            end            
        end
        h=findobj('tag','gbExtraConstraints');
        setappdata(h,'PMFalgIDs',PMFalgIDs);
        flag=get(hcbUpdateMode,'value');
        delete(hmainfig);
        if flag
            commandwindow;
            pause(0.1);
            refineTucker;
        else
            tdalab;
        end
    end   
end