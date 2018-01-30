function PMFActParas=PMFsetOpts(PMFalgList,actModes)
%PMFalgList=[1 3]; for debug
global ScreenWidth ScreenHeight YSpread backGroundColor lFontSize defaultFontName;
MAXPARANUM=10;
currAlg=[];
currIndex=[];
[paraList algs paraTypeList]=PMFalgInit();

algNames={algs(PMFalgList).details};
dParaNames=[];
strParaNames=[];
TFParaNames=[];
hplOpts=[];

PMFActParas=cell(numel(PMFalgList),1);
for i=1:numel(PMFalgList)
    PMFActParas{i}=struct();
end

%% mainWindow
hmainWin=figure('Units','Characters','Resize','on','toolbar','none','menu','none',...
    'tag','PMFmainOpts','Windowstyle','modal','CloseRequestFcn',@optsFinish,...
    'Name','PMF-Options','numbertitle','off','color',backGroundColor);
bkcolor=get(hmainWin,'color');

htextTitle=uicontrol('parent',hmainWin,'Units','normalized','backgroundcolor',bkcolor,'fontsize',lFontSize+4,'fontname',defaultFontName,...
    'string','Options','style','text','fontweight','bold');

hpbOK=uicontrol('parent',hmainWin,'Units','normalized','fontsize',lFontSize,'fontname',defaultFontName,...
    'string','Finish','style','pushbutton','callback',@optsFinish);

hplOpts=uipanel('parent',hmainWin,'Units','normalized',...
    'backgroundcolor',bkcolor,'foregroundcolor','k','borderwidth',1);
htextAlg=uicontrol('parent',hmainWin,'Units','normalized','backgroundcolor',bkcolor,'fontsize',lFontSize,'fontname',defaultFontName,...
    'string',['Algorithm for mode-:' num2str(actModes(1))],'style','text','horizontalalignment','right');

hpmAlg=uicontrol('parent',hmainWin,'Units','normalized','fontsize',lFontSize+2,'fontname',defaultFontName,'foregroundcolor','yellow','backgroundcolor',[ 0    0.4980    0.3294],...
    'string',algNames,'style','popupmenu','callback',@selCurrAlg,'createFcn',@selCurrAlg);

movegui(hmainWin,'center');
waitfor(hmainWin);

    function selCurrAlg(hObject,eventdata,source)
        currIndex=get(hObject,'value');
        currAlg=PMFalgList(currIndex);
        
        set(htextAlg,'string',['Algorithm for mode-',num2str(actModes(currIndex))]);
        
        paraNames=fieldnames(paraList{currAlg});
        inum=min(numel(paraNames),MAXPARANUM);
        paraNames=paraNames(1:inum);

        %% classify the paras
        typeParas=struct2cell(paraTypeList{currAlg});
        i=cellfun(@(x) x==paraType.paraDouble,typeParas);
        dParaNames=paraNames(i);
        numdPara=numel(dParaNames);
        i=cellfun(@(x) x==paraType.paraString,typeParas);
        strParaNames=paraNames(i);
        numstrPara=numel(strParaNames);
        i=cellfun(@(x) x==paraType.paraTF,typeParas);
        TFParaNames=paraNames(i);
        numTFPara=numel(TFParaNames);
        
        lines=ceil(inum-numTFPara/2);

        %% mainWindow
        Width=ScreenWidth/3;
        totalLines=lines+4.5;
        lineHn=1/totalLines;
        lineH=3*YSpread;
        Height=lineH*totalLines;
        nBorderWidth=0.01;
        dnBorderWidth=2*nBorderWidth;
        set(hmainWin,'position',[(ScreenWidth-Width)/2,(ScreenHeight-Height)/2-2*YSpread,Width,Height]);
        set(hplOpts,'parent',hmainWin,'position',[nBorderWidth lineHn 1-dnBorderWidth (lines+1)*lineHn]);
        delete(allchild(hplOpts));
        set(htextAlg,'Position',[nBorderWidth 1-2.2*lineHn 0.3 0.7*lineHn]);
        set(hObject,'Position',[0.32 1-2.4*lineHn 1-dnBorderWidth-0.32 lineHn]);
        set(htextTitle,'Position',[nBorderWidth 1-1.5*lineHn 1-dnBorderWidth 1.1*lineHn]);
        set(hpbOK,'position',[0.4 0.01 0.2 .85*lineHn]);
        
        
        baseline=.95;
        hopt=.95*baseline/lines;
        height=0.7*hopt;
        
        for i=1:numdPara
            uicontrol('parent',hplOpts,'Units','normalized','fontsize',lFontSize,'position',[nBorderWidth,baseline-hopt,0.4,height],...
                'fontname',defaultFontName,'backgroundcolor',bkcolor,'horizontalalignment','right','string',strcat(dParaNames{i},':'),'style','edit','enable','inactive','max',1,'min',0);    
            if isfield(PMFActParas{currIndex},(dParaNames{i}))
                uivalue=PMFActParas{currIndex}.(dParaNames{i});
            else
                uivalue=paraList{currAlg}.(dParaNames{i});
            end
            uicontrol('parent',hplOpts,'Units','normalized','fontsize',lFontSize,'position',[0.42,baseline-hopt,1-dnBorderWidth-0.42,height],'userdata',i,'ButtonDownFcn',@clearField,...
                'fontname',defaultFontName,'horizontalalignment','left','string',uivalue,'style','edit','callback',@getdPara,'createFcn',@getdPara);  
            baseline=baseline-hopt;
        end
        
        for i=1:numstrPara            
            uicontrol('parent',hplOpts,'Units','normalized','fontsize',lFontSize,'position',[nBorderWidth,baseline-hopt,0.4,height],...
                'fontname',defaultFontName,'backgroundcolor',bkcolor,'horizontalalignment','right','string',strcat(strParaNames{i},':'),'style','edit','enable','inactive','max',1,'min',0);  
            
            if isfield(PMFActParas{currIndex},(strParaNames{i}))
                uivalue=PMFActParas{currIndex}.(strParaNames{i});
            else
                uivalue=paraList{currAlg}.(strParaNames{i});
            end
            uicontrol('parent',hplOpts,'Units','normalized','fontsize',lFontSize,'position',[0.42,baseline-hopt,1-dnBorderWidth-0.42,height],'userdata',i,'ButtonDownFcn',@clearField,...
                'fontname',defaultFontName,'horizontalalignment','left','string',uivalue,'style','popupmenu','fontname',defaultFontName,'fontsize',lFontSize+2,'callback',@getstrPara,'createFcn',@getstrPara);  
            baseline=baseline-hopt;
        end
        
        for i=1:numTFPara
            ind=ceil(i/2);
            if isfield(PMFActParas{currIndex},(TFParaNames{i}))
                uivalue=PMFActParas{currIndex}.(TFParaNames{i});
            else
                uivalue=paraList{currAlg}.(TFParaNames{i});
            end
            uicontrol('parent',hplOpts,'Units','normalized','fontsize',lFontSize,'position',[rem(i+1,2)*0.5+dnBorderWidth,baseline-ind*hopt,0.5-dnBorderWidth,height],'userdata',i,...
                'fontname',defaultFontName,'horizontalalignment','left','string',TFParaNames{i},'style','checkbox','fontname',defaultFontName,'fontsize',lFontSize+2,'backgroundcolor',bkcolor,'value',uivalue,'callback',@getTFPara,'createFcn',@getTFPara);  
        end
    end



%% callbacks
    function getdPara(hObject,eventdata,source)
        PMFActParas{currIndex}.(dParaNames{get(hObject,'userdata')})=str2num(get(hObject,'string'));
    end
    function getstrPara(hObject,eventdata,source)
        strpara=get(hObject,'string');
        PMFActParas{currIndex}.(strParaNames{get(hObject,'userdata')})=strtrim(strpara(get(hObject,'value'),:));
    end
    function getTFPara(hObject,eventdata,source)
        PMFActParas{currIndex}.(TFParaNames{get(hObject,'userdata')})=(get(hObject,'value')==1);
    end
    function clearField(hObject,eventdata,source)
        set(hObject,'string','');
    end
    function optsFinish(hObject,eventdata,source)
         PMFActParas=PMFParaListParsing(PMFalgList,PMFActParas,paraList,paraTypeList);
%      save PMFconfig.mat paraList;
        delete(hmainWin);
    end

%% extra functions
    

end

