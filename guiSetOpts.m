function [ opts Main_fig] = guiSetOpts(algs,paraTList,paraList,parasRef)
%% GUI for setting paprmeters
%      algs: structure of algorithms (have fields: name,details
% paraTList: type def
%  paraRefs: default values
global TDALABHOME;
FileDir=horzcat(TDALABHOME,filesep,'userdata',filesep);  %% can be replaced by anyother directory;

if exist('enumeration','file')
    [paraTypes paraTypeNames]=enumeration('paraType');
else
    paraTypes=[paraType.paraDouble paraType.paraMenu paraType.paraTF paraType.paraFile paraType.paraPMF];
    paraTypeNames={'paraDouble' 'paraString' 'paraMenu' 'paraTF' 'paraFile' 'paraPMF'};
end
NPTYPES=numel(paraTypeNames); %% number of types in paraType
hasPAGE2=false;

%% define the pages
maxLs=0;
if iscell(paraList)
    for idx=1:numel(paraList)
        maxLs=max(maxLs,numel(fieldnames(paraList{idx})));
    end
else
    maxLs=numel(fieldnames(paraList));
end
MAX_PER_PAGE=max(5,ceil(maxLs/2));

%% for debug
% algs(1).details='Low-rank NCP';
% algs(1).name='appr3d';
% paraTList{1}=struct('NumOfComp',paraType.paraDouble,'check_5',paraType.paraTF,'MaxIter',paraType.paraDouble,'Tol',paraType.paraDouble,'init',paraType.paraMenu,'initFile',paraType.paraFile,...
%     'PMF',paraType.paraPMF,'check_1',paraType.paraTF,'check_2',paraType.paraTF,'check_3',paraType.paraTF,'check_4',paraType.paraTF);
% paraList{1}=struct('NumOfComp','Number of Components','check_5',false,'MaxIter',100,'Tol',1e-3,'init','random|svds|eigs','initFile','','PMF',[1 2 3],...
%     'check_1',true,'check_2',false,'check_3',true,'check_4',true);
% 
% algs(2).details='Algorithm 2';
% algs(2).name='lraNMF2';
% paraTList{2}=struct('NumOfComp',paraType.paraDouble,'MaxIter',paraType.paraDouble,'Tol',paraType.paraDouble,'init',paraType.paraMenu);
% paraList{2}=struct('NumOfComp','Number of Components','MaxIter',100,'Tol',1e-3,'init','random|svd|eig');



hpage1=[];hpage2=[];hselPage=[];hparas=[];
nAlgs=numel(algs);
currAlg=1;


%% initialization
sortParas;
for n=1:nAlgs
    opts{n}=struct();
end
if ~exist('parasRef','var')
    for n=1:nAlgs
        parasRef{n}=struct();
    end
    copy2opts(paraList);
else
    copy2opts(parasRef);
end


%%
bgColor=[0.8314    0.8157    0.7843];

%% Basic size defs;
Hctrl=3;
Hsplit=.5*Hctrl;

us=get(0,'units');
set(0,'units','characters');
winMain=get(0,'ScreenSize');
WinWidth=winMain(3)/3;
WinHeight=(4.5+MAX_PER_PAGE)*Hctrl+4*Hsplit;
set(0,'units',us);
h1=(2.5*Hctrl+3*Hsplit);

% 

%% main window;
Main_fig=figure('Units','characters',...
    'Resize','on',...
    'toolbar','none',...
    'menu','none',...
    'position',[0 0 WinWidth WinHeight],...
    'Name','Algorithm Options',...
    'numbertitle','off',...
    'color',bgColor,...
    'tag','algOptions','windowstyle','modal',...
    'visible','off');
movegui(Main_fig,'center');


%%Left window
hpnlHelp=uipanel('parent',Main_fig,'units','normalized','backgroundcolor',bgColor,...
    'title','Help','borderwidth',1,...
    'position',[.55 0 .44 1],'visible','off','fontunits','points','fontsize',11,'fontweight','normal');
hpnlPara=uipanel('parent',Main_fig,'units','normalized','backgroundcolor',bgColor,...
    'title','','borderwidth',0,...
    'position',[0 0 1 1]);

%% right
nBorderWidth=.05;
htxtHelp=uicontrol('parent',hpnlHelp,'Units','normalized','position',[0 0.0,1,...
                0.98],'string','','Max',10,'Min',0,...
                'horizontalalignment','left','style','edit','fontunits','points','fontsize',11,'Enable','on');

%% left
hpnlTitle=uipanel('parent',hpnlPara,'units','normalized','backgroundcolor',bgColor,...
    'title','','borderwidth',0,...
    'position',[0 1-h1/WinHeight 1 h1/WinHeight]);
hpage2=uipanel('parent',hpnlPara,'units','normalized','backgroundcolor',bgColor,...
    'title','','borderwidth',0,'highlightcolor','k','BorderType','line','tag','optsPage2','visible','off',...
    'position',[0.025 (2*Hsplit+Hctrl)/WinHeight .95 (MAX_PER_PAGE+1)*Hctrl/WinHeight]);
hpage1=uipanel('parent',hpnlPara,'units','normalized','backgroundcolor',bgColor,...
    'title','','borderwidth',0,'highlightcolor','k','BorderType','line','tag','optsPage1',...
    'position',[0.025 (2*Hsplit+Hctrl)/WinHeight .95 (MAX_PER_PAGE+1)*Hctrl/WinHeight]);
pos=get(hpage1,'position');
hselPage=uicontrol('parent',hpnlPara,'units','normalized',...
    'style','slider','Min',1,'Max',2,'sliderstep',[1 1],'value',2,...
    'position',[pos(3)+pos(1)-4.2/WinWidth (2*Hsplit+1.5*Hctrl+.1)/WinHeight 4/WinWidth ((MAX_PER_PAGE)*Hctrl-.5)/WinHeight],...
    'callback',@selPages);
hpnlButtons=uipanel('parent',hpnlPara,'units','normalized','backgroundcolor',bgColor,...
    'title','','borderwidth',0,...
    'position',[0.01 0.01 0.98 (Hctrl+Hsplit)/WinHeight]);
uicontrol('parent',hpnlPara,'units','normalized','backgroundcolor',bgColor,...
    'style','checkbox','string','Help>>','position',[0.79 (.5*Hsplit+Hctrl)/WinHeight .2 .8*Hctrl/WinHeight] ,...
    'FontUnits','normalized','FontSize',.4,'horizontalalign','right','callback',@showHelp);


%% Title section
uicontrol('parent',hpnlTitle,'Units','normalized','backgroundcolor',bgColor,...
    'style','text','string',' Options ','position',[0. (2*Hsplit+Hctrl)/h1 1 1.2*Hctrl/h1],...
    'FontUnits','normalized','FontSize',.6);
uicontrol('parent',hpnlTitle,'Units','normalized','backgroundcolor',bgColor,...
    'style','text','string','Algorithm:  ','position',[0.05 Hsplit/h1 .25 .8*Hctrl/h1-.02],...
    'FontUnits','normalized','FontSize',.6,'fontweight','normal','horizontalalign','right');
% acolor=[0.7647    0.8392    0.6078];
acolor=[0.9216    0.9451    0.8706];
uicontrol('parent',hpnlTitle,'Units','normalized','backgroundcolor',acolor,...
    'style','popupmenu','string',{algs(:).details},'position',[0.31 Hsplit/h1 .64 .8*Hctrl/h1],...
    'FontUnits','normalized','FontSize',.4,'fontweight','normal','callback',@selAlgs);

%% commands
width=.15;
fontsize=.4;
hpbOK=uicontrol('parent',hpnlButtons,'Units','normalized',...
    'style','pushbutton','string','OK','position',[.05 .15 width .5],...
    'FontUnits','normalized','FontSize',fontsize,'callback',@pbOK);
hpbCancel=uicontrol('parent',hpnlButtons,'Units','normalized',...
    'style','pushbutton','string','Cancel','position',[.23 .15 width .5],...
    'FontUnits','normalized','FontSize',fontsize,'callback',@pbCancel);
hpbDefault=uicontrol('parent',hpnlButtons,'Units','normalized',...
    'style','pushbutton','string','Default','position',[.41 .15 width .5],...
    'FontUnits','normalized','FontSize',fontsize,'callback',@pbDefault);
hpbSave=uicontrol('parent',hpnlButtons,'Units','normalized',...
    'style','pushbutton','string','Save','position',[.59 .15 width .5],...
    'FontUnits','normalized','FontSize',fontsize,'callback',@pbSave);
hpbLoad=uicontrol('parent',hpnlButtons,'Units','normalized',...
    'style','pushbutton','string','Load','position',[.78 .15 width .5],'userdata','NumOfComp',...
    'FontUnits','normalized','FontSize',fontsize,'callback',@pbLoad);

%% Page part
updatePages;
set(Main_fig,'visible','on');
waitfor(Main_fig);



    function selPages(hObject,eventdata,source)
        if get(hObject,'value')==2
            set(hpage1,'visible','on');
            set(allchild(hpage1),'visible','on');
            set(hpage2,'visible','off');
            set(allchild(hpage2),'visible','off');
        else
            set(hpage2,'visible','on');
            set(allchild(hpage2),'visible','on');
            set(hpage1,'visible','off');
            set(allchild(hpage1),'visible','off');
        end        
    end

    function selAlgs(hObject,eventdata,source)
        if currAlg~=get(hObject,'value');
            currAlg=get(hObject,'value');
            updatePages;
        end
    end

    function pbOK(hObject,eventdata,handles)
        defaultParaFix;
        
        %% should be removed
%         for n=1:nAlgs
%             if numel(fieldnames(opts{n}))<1
%                 copy2opts(paraList{n},n);
%             end
%         end
        
%         %% save figs;
%         exportsetupdlg(Main_fig)
%         pause;
   
        delete(Main_fig);
    end

    function pbCancel(hObject,eventdata,handles)
        opts=[];
        delete(Main_fig);
    end

    function pbDefault(hObject,eventdata,handles)
        copy2opts(paraList{currAlg},currAlg);
        updatePages;
    end  

    function pbSave(hObject,eventdata,handles)
        defaultParaFix;
        
%   alg            function_handle              
%   algname        char                         
%   algopts        cell of struct
        algname={algs(:).name};
        alg=cellfun(@str2func,algname,'uni',false);
        algopts=opts;
        
        defname=cellfun(@(x) horzcat(x,'_'),algname,'uni',false);
        defname=cell2mat(defname);
        defname=horzcat(defname,'opts_',date,'.mat');
        
        if numel(algname)==1
            algname=algname{1};
            alg=alg{1};
            algopts=algopts{1};
        end
        [FileName PathName]=uiputfile('*.mat','Save the options as a mat file',horzcat(FileDir,defname));
        if ~isequal(FileName,0)&&~isequal(PathName,0)
            save(horzcat(PathName,FileName),'algname','algopts','alg');
        end
    end

    function pbLoad(hObject,eventdata,handles)
        [FileName PathName]=uigetfile('*.mat','Load settings for algorithms',horzcat(FileDir,'*.mat'));
        if ~isequal(FileName,0)
            fo=load(horzcat(PathName,FileName));
            algNames={algs(:).name};
            if ~iscell(fo.algname)  %% -- to cell
                focell.algname{1}=fo.algname;
                focell.alg{1}=fo.alg;
                focell.algopts{1}=fo.algopts;
                fo=focell;
                clear focell;
            end
            for cidx=1:numel(fo.algname)
                pos=find(strcmpi(algNames,fo.algname{cidx})==1);
                if ~isempty(pos)
                    for iag=1:numel(pos)
                        copy2opts(fo.algopts{cidx},iag);
                    end
                end
            end
        end
        updatePages;
    end

    function sortParas
        %% sort paraTList,paraList
        paraTypeOrder=[paraType.paraDouble paraType.paraString paraType.paraMenu paraType.paraFile paraType.paraPMF paraType.paraTF];
        for na=1:nAlgs
            paraNames=fieldnames(paraTList{na});
            numTypes=struct2cell(paraTList{na});
            L=1:numel(numTypes);
            ord=[];
            for cidx=1:size(paraTypeOrder,2)
                flag=(cellfun(@(x) x==paraTypeOrder(cidx),numTypes));
                ord(end+1:end+sum(flag))=L(flag);
            end
            paraNames=paraNames(ord);
            paraTList{na}=orderfields(paraTList{na},paraNames);
            paraList{na}=orderfields(paraList{na},paraNames);
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%% UPDAGE PAGES 
    function updatePages
        set(hpage2,'visible','off');
        set(hpage1,'visible','off');
        if ~isempty(hparas)
            delete(get(hpage1,'child'));
            delete(get(hpage2,'child'));
            hparas=[];
        end
        hasPAGE2=false;
        
        paraNames=fieldnames(paraTList{currAlg});
        numTypes=struct2cell(paraTList{currAlg});
        NPTypes=numel(paraNames);
        nPTypes=zeros(NPTYPES,1);
        for i=1:NPTYPES
            flag=cellfun(@(x) x==paraType.(paraTypeNames{i}),numTypes);
            nPTypes(i)=sum(flag==1);
        end
        PLIndex=1;  % Line index
        hPIndex=1;  % Para index
        
        Hruler=1/(MAX_PER_PAGE+1);
        Hctrl=Hruler*.65;
        

%% 1.paraType.paraDouble
        for i=1:nPTypes(int8(paraType.paraDouble))
            lino=rem(PLIndex-1,MAX_PER_PAGE)+1.5;
            if PLIndex<=MAX_PER_PAGE
                currPage=hpage1;
            else
                currPage=hpage2;
                hasPAGE2=true;
            end
            
            
            %% uicontrols
            uicontrol('parent',currPage,'units','normalized','position',[0.01 1-lino*Hruler 0.3 Hctrl],'backgroundcolor',bgColor,...
                'style','text','string',horzcat(paraNames{hPIndex},': '),'horizontalalignment','right',...
                'fontunits','normalized','fontsize',.5,'fontweight','normal');
            if isfield(opts{currAlg},paraNames{hPIndex})
                strValue=num2str(opts{currAlg}.(paraNames{hPIndex}));
            else
                strValue=num2str(paraList{currAlg}.(paraNames{hPIndex}));
            end
            hparas(hPIndex)=uicontrol('parent',currPage,'units','normalized','position',[0.35 1.025-lino*Hruler 0.5 Hctrl],...
                'style','edit','string',strValue,'horizontalalignment','left','ButtonDownFcn',@clearField,...
                'fontunits','normalized','fontsize',.5,'userdata',paraNames{hPIndex},'callback',@setparaDouble);
            hPIndex=hPIndex+1;
            PLIndex=PLIndex+1;
        end
        
        
%% 2.paraType.paraString
        for i=1:nPTypes(int8(paraType.paraString))
            lino=rem(PLIndex-1,MAX_PER_PAGE)+1.5;
            if PLIndex<=MAX_PER_PAGE
                currPage=hpage1;
            else
                currPage=hpage2;
                hasPAGE2=true;
            end
            
            %% uicontrols
            uicontrol('parent',currPage,'units','normalized','position',[0.01 1-lino*Hruler 0.3 Hctrl],'backgroundcolor',bgColor,...
                'style','text','string',horzcat(paraNames{hPIndex},': '),'horizontalalignment','right',...
                'fontunits','normalized','fontsize',.5,'fontweight','normal');
            
            hparas(hPIndex)=uicontrol('parent',currPage,'units','normalized','position',[0.35 1.025-lino*Hruler 0.5 Hctrl],...
                'style','edit','string',paraList{currAlg}.(paraNames{hPIndex}),'horizontalalignment','left','ButtonDownFcn',@clearField,...
                'fontunits','normalized','fontsize',.5,'value',1,'userdata',paraNames{hPIndex},'callback',@setparaString);
            if isfield(opts{currAlg},paraNames{hPIndex})
                set(hparas(hPIndex),'string',opts{currAlg}.(paraNames{hPIndex}));
            end
            hPIndex=hPIndex+1;
            PLIndex=PLIndex+1;
        end
        
%% 3.paraType.paraMenu
        for i=1:nPTypes(int8(paraType.paraMenu))
            lino=rem(PLIndex-1,MAX_PER_PAGE)+1.5;
            if PLIndex<=MAX_PER_PAGE
                currPage=hpage1;
            else
                currPage=hpage2;
                hasPAGE2=true;
            end
            
            %% uicontrols
            uicontrol('parent',currPage,'units','normalized','position',[0.01 1-lino*Hruler 0.3 Hctrl],'backgroundcolor',bgColor,...
                'style','text','string',horzcat(paraNames{hPIndex},': '),'horizontalalignment','right',...
                'fontunits','normalized','fontsize',.5,'fontweight','normal');
            
            parastrs=parastrsplit(paraList{currAlg}.(paraNames{hPIndex}));
            hparas(hPIndex)=uicontrol('parent',currPage,'units','normalized','position',[0.35 1.025-lino*Hruler 0.5 Hctrl],...
                'style','popupmenu','string',parastrs,'horizontalalignment','left',...
                'fontunits','normalized','fontsize',.5,'value',1,'userdata',paraNames{hPIndex},'callback',@setparaMenu);
            if isfield(opts{currAlg},paraNames{hPIndex})  
                selValue=find(strcmpi(get(hparas(hPIndex),'string'),opts{currAlg}.(paraNames{hPIndex}))==1);
                set(hparas(hPIndex),'value',selValue);
            end

            hPIndex=hPIndex+1;
            PLIndex=PLIndex+1;
        end
     
%% 4.paraType.paraFile
        for i=1:nPTypes(int8(paraType.paraFile))
            lino=rem(PLIndex-1,MAX_PER_PAGE)+1.5;
            if PLIndex<=MAX_PER_PAGE
                currPage=hpage1;
            else
                currPage=hpage2;
                hasPAGE2=true;
            end
            
            %% uicontrols
            uicontrol('parent',currPage,'units','normalized','position',[0.01 1-lino*Hruler 0.3 Hctrl],'backgroundcolor',bgColor,...
                'style','text','string',horzcat(paraNames{hPIndex},': '),'horizontalalignment','right',...
                'fontunits','normalized','fontsize',.5,'fontweight','normal');
            if isfield(opts{currAlg},paraNames{hPIndex})
                opts{currAlg}.(paraNames{hPIndex})=reshape(opts{currAlg}.(paraNames{hPIndex}),1,[]);
                strValue=opts{currAlg}.(paraNames{hPIndex});
            else
                strValue='Click to load a file.';
            end
            hparas(hPIndex)=uicontrol('parent',currPage,'units','normalized','position',[0.35 1.025-lino*Hruler 0.5 Hctrl],...
                'style','edit','string',strValue,'horizontalalignment','left','ButtonDownFcn',@setparaFile,...
                'fontunits','normalized','fontsize',.5,'enable','inactive','userdata',paraNames{hPIndex});
            hPIndex=hPIndex+1;
            PLIndex=PLIndex+1;
        end
      
%% 5.paraType.paraPMF
        for i=1:nPTypes(int8(paraType.paraPMF))
            lino=rem(PLIndex-1,MAX_PER_PAGE)+1.5;
            if PLIndex<=MAX_PER_PAGE
                currPage=hpage1;
            else
                currPage=hpage2;
                hasPAGE2=true;
            end
            
            %% uicontrols
            uicontrol('parent',currPage,'units','normalized','position',[0.01 1-lino*Hruler 0.3 Hctrl],'backgroundcolor',bgColor,...
                'style','text','string',horzcat(paraNames{hPIndex},': '),'horizontalalignment','right',...
                'fontunits','normalized','fontsize',.5,'fontweight','normal');
            if isfield(opts{currAlg},paraNames{hPIndex})
                strValue=num2str(opts{currAlg}.(paraNames{hPIndex}));
            else
                strValue='Click to select algorithms.';
            end
            hparas(hPIndex)=uicontrol('parent',currPage,'units','normalized','position',[0.35 1.025-lino*Hruler 0.5 Hctrl],...
                'style','edit','string',strValue,'horizontalalignment','left','ButtonDownFcn',@setPMFAlgs,...
                'fontunits','normalized','fontsize',.5,'enable','inactive','userdata',paraNames{hPIndex});
            hPIndex=hPIndex+1;
            PLIndex=PLIndex+1;
        end  
        
         
%% LAST.paraType.paraTF
        for i=1:nPTypes(int8(paraType.paraTF))
            lino=rem(PLIndex-1,MAX_PER_PAGE)+1.5;
            if PLIndex<=MAX_PER_PAGE
                currPage=hpage1;
            else
                currPage=hpage2;
                hasPAGE2=true;
            end
            
            if isfield(opts{currAlg},paraNames{hPIndex})
                tfValue=opts{currAlg}.(paraNames{hPIndex});
            else
                tfValue=paraList{currAlg}.(paraNames{hPIndex});
            end
            
            %% uicontrols
            hparas(hPIndex)=uicontrol('parent',currPage,'units','normalized','position',[.1+rem(i-1,2)*.5 1.025-lino*Hruler 0.35 Hctrl],...
                'style','checkbox','string',paraNames{hPIndex},'horizontalalignment','left','backgroundcolor',bgColor,'callback',@setparaTF,...
                'fontunits','normalized','fontsize',.5,'enable','on','fontweight','normal','value',tfValue,'userdata',paraNames{hPIndex},...
                'Max',true,'Min',false);
            hPIndex=hPIndex+1;
            PLIndex=PLIndex+rem(i-1,2);
        end
        
        if hasPAGE2
            set(hselPage,'value',2,'visible','on');
        else
            set(hselPage,'value',2,'visible','off');
        end
        set(hpage1,'visible','on');
        set(hpage2,'visible','off');
 
    end %% end of updatePages       

%% update parameters ----  copy curropts to the output opts
    function copy2opts(curropts,indices)
        if nargin==1
            indices=1:nAlgs;
        end
        if numel(indices)>1
            for na=indices
                names=fieldnames(curropts{na});
                for cidx=1:numel(names)
                    switch int8(paraTList{na}.(names{cidx}))
                        case {int8(paraType.paraDouble),int8(paraType.paraPMF)}
                            if isnumeric(curropts{na}.(names{cidx}))
                                opts{na}.(names{cidx})=reshape(curropts{na}.(names{cidx}),1,[]);
                            end
                        case {int8(paraType.paraString),int8(paraType.paraFile)}
                            if ischar(curropts{na}.(names{cidx}))
                                opts{na}.(names{cidx})=curropts{na}.(names{cidx});
                            end
                        case int8(paraType.paraMenu)
                            selstrs=parastrsplit(paraList{na}.(names{cidx}));
                            pos=find(strcmpi(selstrs,curropts{na}.(names{cidx}))==1);
                            if ~isempty(pos)
                                opts{na}.(names{cidx})=selstrs{pos};
                            else
                                opts{na}.(names{cidx})=selstrs{1};
                            end
                        case int8(paraType.paraTF)
                                opts{na}.(names{cidx})=logical(curropts{na}.(names{cidx}));
                        otherwise
                            errordlg(horzcat('Unknown data type of parameters [',names{cidx},']'),'Data type error');
                    end % switch end
                end % cidx-names
            end % end of algs
        else % indices==1
            na=indices;
            if iscell(curropts)
                curropts=curropts{1};
            end
            names=fieldnames(curropts);
            for cidx=1:numel(names)
                switch int8(paraTList{na}.(names{cidx}))
                    case {int8(paraType.paraDouble),int8(paraType.paraPMF)}
                        if isnumeric(curropts.(names{cidx}))
                            opts{na}.(names{cidx})=reshape(curropts.(names{cidx}),1,[]);
                        end
                        case {int8(paraType.paraFile),int8(paraType.paraString)}
                            if ischar(curropts.(names{cidx}))
                                opts{na}.(names{cidx})=curropts.(names{cidx});
                            end
                        case int8(paraType.paraMenu)
                            selstrs=parastrsplit(paraList{na}.(names{cidx}));
                            pos=find(strcmpi(selstrs,curropts.(names{cidx}))==1);
                            if ~isempty(pos)
                                opts{na}.(names{cidx})=selstrs{pos};
                            else
                                opts{na}.(names{cidx})=selstrs{1};
                            end
                    case int8(paraType.paraTF)
                            opts{na}.(names{cidx})=logical(curropts.(names{cidx}));
                    otherwise
                        errordlg(horzcat('Unknown data type of parameters [',names{cidx},']'),'Data type error');
                end % switch end
            end % cidx-names
        end
    end

%% defaultParaFix  -- check and update the output opts by using the default values
    function defaultParaFix
        for cidx=1:nAlgs
            defaultNames=fieldnames(paraList{cidx});
            for i=1:numel(defaultNames)
                if ~isfield(opts{cidx},defaultNames{i}) %% trying to fix this parameter
                    switch int8(paraTList{cidx}.(defaultNames{i}))
                        case {int8(paraType.paraDouble),int8(paraType.paraPMF)}
                            value=(paraList{cidx}.(defaultNames{i}));
                            if isnumeric(value)
                                opts{cidx}.(defaultNames{i})=value;
                            end
                        case {int8(paraType.paraMenu)}
                            strs=parastrsplit(paraList{cidx}.(defaultNames{i}));
                            opts{cidx}.(defaultNames{i})=strs{1};
                        case {int8(paraType.paraFile),int8(paraType.paraString)}
                            opts{cidx}.(defaultNames{i})=(paraList{cidx}.(defaultNames{i}));
                        case {int8(paraType.paraTF)}
                            opts{cidx}.(defaultNames{i})=logical(paraList{cidx}.(defaultNames{i}));
                    end %% switch
                else %% check the status
                    switch int8(paraTList{cidx}.(defaultNames{i}))
                        case {int8(paraType.paraDouble),int8(paraType.paraPMF)}
                            opts{cidx}.(defaultNames{i})=double(opts{cidx}.(defaultNames{i}));
                        case {int8(paraType.paraMenu)}
                            if isnumeric(opts{cidx}.(defaultNames{i}))
                                opts{cidx}.(defaultNames{i})=num2str(opts{cidx}.(defaultNames{i}));
                            else
                                opts{cidx}.(defaultNames{i})=char(opts{cidx}.(defaultNames{i}));
                            end
                        case {int8(paraType.paraFile),int8(paraType.paraString)}
                            if isempty(opts{cidx}.(defaultNames{i}))
                                opts{cidx}.(defaultNames{i})='NotSetYet';
                            end
                        case {int8(paraType.paraTF)}
                            opts{cidx}.(defaultNames{i})=logical(opts{cidx}.(defaultNames{i}));
                    end %% switch                    
                end %% not set yet
            end %% check all the field            
            
        end %% for all algorithms
    end


%% supporting functions
    function c=parastrsplit(s) % s is split by '|'. empty elements will be removed from c;
        c=regexp(s,'\s*\|\s*','split');
        c=strtrim(c);
        c=c(~cellfun('isempty', c));
    end

    function clearField(hObject,eventdata,handles)
        set(hObject,'string','');
    end

    function setPMFAlgs(hObject,eventdata,handles)
        idx=setPMFAlgIndices;
        idx=reshape(idx,1,[]);
        if ~isempty(idx)
            set(hObject,'string',num2str(idx));
            opts{currAlg}.(get(hObject,'userdata'))=idx;
        end
    end

    function showHelp(hObject,eventdata,handles)
        if get(hObject,'value')  %% show help
            helpInfor=help(algs(currAlg).name);
            helpInfor=regexp(helpInfor, '\n', 'split')';
            if isempty(helpInfor)
                helpInfor='Help information unavailable.';
            end
            set(htxtHelp,'string',helpInfor);
            
            pos=get(Main_fig,'position');
            pos(3)=pos(3)*2;
            set(Main_fig,'position',pos);
            set(hpnlPara,'position',[0 0 .54 1]);
            set(hpnlHelp,'visible','on','title',horzcat('[Help for the ',algs(currAlg).details,']'));
            movegui(Main_fig,'center');
        else
            pos=get(Main_fig,'position');
            pos(3)=pos(3)./2;
            set(Main_fig,'position',pos);
            set(hpnlPara,'position',[0 0 1 1]);
            set(hpnlHelp,'visible','off');  
            movegui(Main_fig,'center');
            set(htxtHelp,'string','');
        end
    end

%% for paras callback
    function setparaDouble(hObject,eventdata,handles)
        opts{currAlg}.(get(hObject,'userdata'))=str2num(get(hObject,'string'));
    end
    function setparaString(hObject,eventdata,handles)
        opts{currAlg}.(get(hObject,'userdata'))=get(hObject,'string');
    end
    function setparaMenu(hObject,eventdata,handles)
        strs=(get(hObject,'string'));
        opts{currAlg}.(get(hObject,'userdata'))=strs{get(hObject,'value')};
    end
    function setparaFile(hObject,eventdata,handles)
        [FileName FilePath]=uigetfile('*.*','Select a file',horzcat(FileDir,'*.mat'));
        if ~isequal(FileName,0)
            FileName=fullfile(FilePath,FileName);
            set(hObject,'string',FileName);
            opts{currAlg}.(get(hObject,'userdata'))=FileName;
        end
    end
    function setparaTF(hObject,eventdata,handles)
        opts{currAlg}.(get(hObject,'userdata'))=logical(get(hObject,'value'));
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF THE FUNCTION
end

