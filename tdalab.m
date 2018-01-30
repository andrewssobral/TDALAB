function tdalab(action)
%% Tensor decomposition and analysis LAB
%%  by Guoxu Zhou
%%  ver: 2011.4.1


%% Only one instance is permitted
if nargin==0
    action='show';
end
h=allchild(0);
oo = findall(h, 'Tag', 'TDALAB' );
if ~isempty(oo);
    if strcmpi(action,'hide')
        set(oo,'visible','off','HandleVisibility','off');
    else
        set(oo,'visible','on','HandleVisibility','callback');
        figure(oo);
    end
    return;
end

%% general initialization
clear all;
clc;


%% global variables
global TDALABPATH TEMPFILE  TDALABHOME;
global haxesTensorInfor hmainTDALAB;
global ScreenWidth ScreenHeight;
global nFontSize defaultFontName sFontSize lFontSize;
global mainWidth mainHeight;
global backGroundColor pbBackGroundColor;
global algs tdalabStatus Y Ycap NumOfComp NumOfMode;


%% initialization
TDALABconfig;
%% Tensor toolbox is required.
if any([~exist('tensor','file'), ~exist('ktensor','file'), ~exist('ttensor','file')])
    disp('Tensor Toolbox ver2.4 or above is required.');
    fprintf('Download <a href = "http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html">Tensor toolbox Ver2.5</a> and/or check your path settings.\n');
    rmpath(TDALABPATH);
    return;
end
nBorderWidth=0.01;dnBorderWidth=2*nBorderWidth;

%% load algorithms information
TDALABinitialization;

%% constant variables
inputFormat=defstr('input');
TDModelstr=defstr('TDModel');
ModeConstr=defstr('TuckerConstr');

%% defin mainWindow
mainHeight=min(65,ScreenHeight-5);mainWidth=min(130,ScreenWidth);
hmainTDALAB=figure('Units','Characters','Resize','on','toolbar','none','menu','none',...
    'OuterPosition',[0,0,mainWidth,mainHeight],...
    'tag','TDALAB','Name',defstr('TDALAB'),'numbertitle','off','color',backGroundColor,'visible','off','CloseRequestFcn',@TDALABexit);
pos=get(hmainTDALAB,'Position');mainWidth=pos(3);mainHeight=pos(4);
colormap('default');
bkcolor=get(hmainTDALAB,'Color');

%% menu
hmenuFile=uimenu('label','&File','tag','menuFile');
huimenuApp=uimenu('label','&Applications','tag','menuApplications');
% hmenuDemo=uimenu('label','&DEMO','tag','menuDemo');
hmenuHelp=uimenu('label','&Help','tag','menuHelp');
uimenu(hmenuFile,'label','Load','tag','menuFile_Open','callback','loadTensor;');
uimenu(hmenuFile,'label','Save Results','tag','menuFile_Export','callback',@saveYcapOnly);
uimenu(hmenuFile,'label','Save Workspace','tag','menuFile_SaveAs','Separator','on','callback','saveWorkSpace;');
uimenu(hmenuFile,'label','Load Workspace','tag','menuFile_WSLoad', 'Separator','off','callback','LoadWorkSpace;');
uimenu(hmenuFile,'label','Export Setup...','tag','menuFile_Export','Separator','on','callback',@FigureExport);
uimenu(hmenuFile,'label','Exit','tag','menuFile_Exit','Separator','on','callback',@TDALABexit);
uimenu(huimenuApp,'label','Tucker discriminant analysis','tag','menuApp_TensorDA','callback','callTuckerDA;');
uimenu(huimenuApp,'label','Clustering analysis','tag','menuApp_cluster','callback','runclustering;');


hmenuHelp_update=uimenu(hmenuHelp,'label','Update...','callback','TDALAB_Update');
hmenuHelp_about=uimenu(hmenuHelp,'label','About','callback','TDALAB_About');

%% logo
hLogo=uicontrol('parent',hmainTDALAB,'Units','Normalized','position',[0. 0.9 1 0.09],'BackgroundColor',bkcolor,'style','text',...
    'string','T  D  A  L  A  B','FontUnits','Normalized','FontSize',nFontSize,'HandleVisibility','off','backgroundcolor',bkcolor);
uicontrol('parent',hmainTDALAB,'Units','Normalized','position',[0 0.85 1 0.02],'style','text','string',...
    '[ L a b o r a t o r y     f o r    T e n s o r    D e c o m p o s i t i o n    &    A n a l y s i s ]',...
    'FontUnits','normalized','FontSize',nFontSize,'BackgroundColor',bkcolor,'FontName',defaultFontName);


%% Load tensor
hpanelProb=uipanel('parent',hmainTDALAB,'Units','normalized','BorderType','line','position',...
    [nBorderWidth 0.57 0.55-nBorderWidth 0.25],'backgroundcolor',bkcolor,'title','',...
    'FontUnits','points','fontsize',sFontSize,'foregroundcolor','k','borderwidth',1,'highlightcolor','k');

uicontrol('parent',hpanelProb,'Units','normalized','style','text','fontunits','points','fontsize',lFontSize,'string',...
    'Select variable:','position',[nBorderWidth,0.85,0.35,0.08],'horizontalalign','left','fontname',defaultFontName,...
    'backgroundcolor',bkcolor,'enable','off','tag','txtSelVar');
uicontrol('parent',hpanelProb,'Units','normalized','style','popupmenu','fontsize',lFontSize,'fontunits','points','string',' ',...
    'position',[0.5,0.86,0.5-nBorderWidth,0.08],'horizontalalign','left',...
    'fontname',defaultFontName,'enable','on','tag','ppSelVar','callback','selectVars;');
%
uicontrol('parent',hpanelProb,'Units','normalized','style','text','fontunits','points','fontsize',lFontSize,'string','Mode:',...
    'position',[nBorderWidth,0.66,0.15,0.08],'horizontalalign','left','fontname',defaultFontName,...
    'backgroundcolor',bkcolor,'enable','off','tag','txtMode');
uicontrol('parent',hpanelProb,'Units','Normalized','style','popupmenu','fontsize',lFontSize+1,'string',' ',...
    'position',[nBorderWidth+0.15,0.685,0.3,0.08],'horizontalalign','left',...
    'backgroundcolor','white','fontname',defaultFontName,'enable','on','tag','ppMode','enable','off');
uicontrol('parent',hpanelProb,'Units','normalized','style','pushbutton','fontunits','points','fontsize',lFontSize,'string','View & Modify',...
    'position',[0.5,0.62,0.5-nBorderWidth,0.14],'horizontalalign','left','backgroundcolor',pbBackGroundColor,...
    'fontname',defaultFontName,'enable','off','tag','pbModi','callback','updateMode;');
%
haxesTensorInfor=axes('Units','normalized','position',[dnBorderWidth,0.58,0.54-nBorderWidth,0.12],...
    'Color','none','visible','off');
uicontrol('parent',hpanelProb,'Units','normalized','style','text','backgroundcolor',bkcolor,...
    'position',[0.02,0.02,0.96,0.48],'tag','txtNoTensorInfor','string',{'Please load a tensor first.','[Menu-File-Open]'},...
    'horizontalAlign','left','Fontsize',lFontSize,'FontName',defaultFontName,'Max',10,'Min',1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tucker core / noise
h=0.125*0.8-nBorderWidth/2;
hpanelInputFormat=uipanel('parent',hmainTDALAB,'Units','normalized','BorderType','line','borderWidth',1,...
    'position',[0.55+nBorderWidth 0.57+h+nBorderWidth 0.45-dnBorderWidth 0.25-h-nBorderWidth],'backgroundcolor',bkcolor,'highlightcolor','k');
%% input
uicontrol('parent',hpanelInputFormat,'Units','normalized','position',[nBorderWidth,0.59,0.4-dnBorderWidth,0.25],'style','text','string','Input format:','fontunits','points','horizontalalign','left','backgroundcolor',bkcolor,...
    'tag','txtInputFormat','fontsize',lFontSize,'fontname',defaultFontName,'enable','off');
uicontrol('parent',hpanelInputFormat,'Units','normalized','position',[0.4,0.635,0.6-dnBorderWidth,0.25],'style','popupmenu','string',defstr('input'),'enable','off',...
    'TooltipString','ktensor/ttensor','tag','pmInputFormat','fontsize',lFontSize,'fontname',defaultFontName,'callback',@selInputFormat);
uicontrol('parent',hpanelInputFormat,'Units','normalized','style','text','fontsize',lFontSize,'string','Core tensor:',...
    'position',[nBorderWidth 0.285 0.4-dnBorderWidth 0.25],'horizontalalign','left','fontname',defaultFontName,...
    'backgroundcolor',bkcolor,'enable','off','tag','txtCoreTensor');
GGenMode=defstr('gencore');
uicontrol('parent',hpanelInputFormat,'Units','normalized','style','popupmenu','fontsize',lFontSize,'string',GGenMode,...
    'position',[0.4 0.3 0.6-dnBorderWidth 0.25],'horizontalalign','left','fontname',defaultFontName,...
    'enable','off','tag','pmG','callback',@genCoreTensor);
%% full tensor
hckToTensor=uicontrol('parent',hpanelInputFormat,'Units','normalized','style','checkbox','fontsize',lFontSize,'string','Full tensor',...
    'position',[0.6,0.01,0.4-nBorderWidth,0.25],'HorizontalAlignmen','right','fontname',defaultFontName ,...
    'enable','on','tag','ckToTensor','backgroundColor',bkcolor,'callback',@isFullTensor,'enable','off','TooltipString',...
    'Using tensor as input instead of decomposed form.');


hpanelNoise=uipanel('parent',hmainTDALAB,'Units','normalized','BorderType','line','borderWidth',1,...
    'position',[0.55+nBorderWidth 0.57 0.45-dnBorderWidth h],'backgroundcolor',bkcolor);
set(hpanelNoise,'highlightcolor','k');
htxtNoise=uicontrol('parent',hpanelNoise,'Units','normalized','string','Add noise',...
    'position',[nBorderWidth 0.5 1-dnBorderWidth 0.4],'backgroundcolor',bkcolor,...
    'fontname',defaultFontName,'fontsize',lFontSize,'style','text','tag','txtNoise','enable','off');
NoiseStr={'no noise';'Uniform';'Gaussian';'Abs. Gaussian'};
uicontrol('parent',hpanelNoise,'Units','normalized','string',NoiseStr,'enable','off',...
    'position',[nBorderWidth 0.1 1-dnBorderWidth 0.4],'tag','pmNoise',...
    'fontname',defaultFontName,'fontsize',lFontSize,'style','popupmenu','tag','pmNoise','callback',@setNoise);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select Algorithm
hpanelAlg=uipanel('parent',hmainTDALAB,'Units','normalized','BorderType','line','position',...
    [nBorderWidth 0.33 1-dnBorderWidth 0.24-nBorderWidth],'backgroundcolor',bkcolor,...
    'fontsize',lFontSize,'fontname',defaultFontName,'borderwidth',1,'highlightcolor','k');

modelTipstr=sprintf('   CP: CANDECOMP/PARAFAC model \nTucker: Tucker model \n  PMF: Penalized Matrix Factorization (BSS, ICA, NMF, ...) \n  BCD: Block Component/Term Decomposition');
uicontrol('parent',hpanelAlg,'Units','normalized','fontsize',lFontSize,'string',TDModelstr,'TooltipString',modelTipstr,...
    'position',[0.18 0.8 0.2-nBorderWidth 0.15],'fontname',defaultFontName,...
    'style','popupmenu','visible','on','callback',@selTDModel,'tag','pmTDModel');

uicontrol('parent',hpanelAlg,'Units','normalized','string','TD model:','Fontname',defaultFontName,...
    'style','text','position',[nBorderWidth 0.77 0.15-nBorderWidth 0.15],'fontsize',lFontSize,'horizontalalignment','left',...
    'backgroundcolor',bkcolor,'tag','txtTDModel','enable','on');

% NonNegativity
uicontrol('parent',hpanelAlg,'Units','normalized','string','Nonnegative','Fontname',defaultFontName,...
    'style','checkbox','position',[0.5 0.8 0.2 0.15],'fontsize',lFontSize,'value',tdalabStatus.nonnegativity,...
    'backgroundcolor',bkcolor,'max',true,'min',false,'tag','cbNN','enable','on','callback',...
    @setNonegative);

%% ExtraConstraints for tucker
uicontrol('parent',hpanelAlg,'Units','normalized','string','[Tucker]Feature extraction:','Fontname',defaultFontName,...
    'style','text','position',[nBorderWidth 0.55 0.28-nBorderWidth 0.12],'fontsize',lFontSize,'horizontalalignment','left',...
    'backgroundcolor',bkcolor,'enable','on');
hgbExtraConstraints=uibuttongroup('parent',hpanelAlg,'Units','normalized','fontsize',lFontSize,...
    'position',[0.3 0.55 0.7-dnBorderWidth 0.15],'backgroundcolor',bkcolor,...
    'bordertype','line','borderwidth',0,'visible','on','SelectionChangeFcn',@setExtraCons,'tag','gbExtraConstraints');
% width=(1-dnBorderWidth)/numel(ModeConstr);
strlen=max(cellfun(@length,ModeConstr),8);
width=strlen/sum(strlen);
for i=1:numel(ModeConstr)
    uicontrol('parent',hgbExtraConstraints,'Units','normalized','fontsize',lFontSize,...
        'position',[nBorderWidth+sum(width(1:i-1)) 0 width(i) 1],'fontname',defaultFontName,'backgroundcolor',bkcolor,...
        'style','radiobutton','tag','rbExtraConstraints','string',ModeConstr{i},'userdata',i,'enable','off');
end
set(hgbExtraConstraints,'SelectedObject',[]);

htextAlgs=uicontrol('parent',hpanelAlg,'Units','normalized','string','Algorithm List:',...
    'style','text','position',[nBorderWidth 0.26 0.15-nBorderWidth 0.15],'fontsize',lFontSize,...
    'backgroundcolor',bkcolor,'horizontalalign','left','fontname',defaultFontName);
%% define name list for popupmenu
algList={algs(:).details};
uicontrol('parent',hpanelAlg,'Units','normalized','string','Algorithm List:','fontname',defaultFontName,...
    'style','popupmenu','position',[0.16 0.3 0.68 0.15],'fontsize',lFontSize,...
    'tag','pmAlgList','string',algList,'callback',@setAlg);
uicontrol('parent',hpanelAlg,'Units','normalized','string','Show all','fontname',defaultFontName,...
    'style','checkbox','position',[0.86 0.3 0.14-nBorderWidth 0.15],'fontsize',lFontSize,...
    'tag','cbShowAllAlgs','callback','updateAlgList;','backgroundcolor',bkcolor,'value',false);

updateAlgList;
hpbOpts=uicontrol('parent',hpanelAlg,'Units','normalized','string','Advanced Options','fontname',defaultFontName,...
    'style','pushbutton','position',[0.75 4*nBorderWidth 0.25-nBorderWidth 0.2],'fontsize',lFontSize,...
    'tag','pbOpts','callback','setOpts;','enable','on','backgroundcolor',pbBackGroundColor);

uicontrol('parent',hpanelAlg,'Units','normalized','string','NumOfComp','fontname',defaultFontName,...
    'style','checkbox','position',[nBorderWidth 4*nBorderWidth 0.25-nBorderWidth 0.22],'fontsize',lFontSize,...
    'tag','chkmoderank','callback',@setNumOfComp,'enable','on','backgroundcolor',backGroundColor,'max',true,'min',false);

%% Run The Selected Algorithm
hpbRunAlg=uicontrol('parent',hmainTDALAB,'Units','normalized','string','Run Now!','backgroundcolor',pbBackGroundColor,...
    'style','pushbutton','position',[0.31 0.265 0.24 0.05],'fontsize',lFontSize,'enable','off',...
    'tag','pbRunAlg','callback','runTDalg');

%% Output analysis
hpanelOutput=uipanel('parent',hmainTDALAB,'Units','normalized','BorderType','line','position',...
    [nBorderWidth 0.08 0.55-nBorderWidth 0.17],'backgroundcolor',bkcolor,'title','',...
    'FontUnits','points','fontsize',sFontSize,'foregroundcolor','k','tag','pnlOutputAnalysis');
set(hpanelOutput,'borderwidth',1,'highlightcolor','k');
htxtOutput=uicontrol('parent',hpanelOutput,'Units','normalized','style','text','backgroundcolor',bkcolor,...
    'position',[nBorderWidth 0.7 1-dnBorderWidth 0.2],'string','Output analysis','FontSize',lFontSize,...
    'FontName',defaultFontName,'enable','off');
w=0.38-dnBorderWidth;
hpbVisual=uicontrol('parent',hpanelOutput,'Units','normalized','style','pushbutton','backgroundcolor',pbBackGroundColor,...
    'position',[0.1 0.4 w 0.25],'string','Visualization','FontSize',lFontSize,'callback','visualize;',...
    'FontName',defaultFontName,'enable','off');
hpbCloseFigs=uicontrol('parent',hpanelOutput,'Units','normalized','style','pushbutton','backgroundcolor',pbBackGroundColor,...
    'position',[0.1 0.05 w 0.25],'string','Close figures','FontSize',lFontSize,...
    'FontName',defaultFontName,'callback',@closeFigs,'enable','off');
hpbCCW=uicontrol('parent',hpanelOutput,'Units','normalized','style','pushbutton','backgroundcolor',pbBackGroundColor,...
    'position',[0.9-w 0.4 w 0.25],'string','Command win>>','FontSize',lFontSize,...
    'FontName',defaultFontName,'callback',@gotoCCW,'enable','on','tag','pbgotoCCwin');
uicontrol('parent',hpanelOutput,'Units','normalized','style','pushbutton','backgroundcolor',pbBackGroundColor,...
    'position',[0.9-w 0.05 w 0.25],'string','Save','FontSize',lFontSize,...
    'FontName',defaultFontName,'callback',@saveYcapOnly,'enable','off');
hpanelAdvComp=uipanel('parent',hmainTDALAB,'Units','normalized','BorderType','line','position',...
    [0.55+nBorderWidth 0.08 0.45-dnBorderWidth 0.17],'backgroundcolor',bkcolor,'title','','tag','pnlAdvComp',...
    'FontUnits','points','fontsize',sFontSize,'foregroundcolor','k','borderwidth',1,'highlightcolor','k');
uicontrol('parent',hpanelAdvComp,'Units','normalized','style','text','backgroundcolor',bkcolor,...
    'position',[nBorderWidth 0.7 1-dnBorderWidth 0.2],'string','Algorithm analysis','FontSize',lFontSize,...
    'FontName',defaultFontName,'enable','off','tag','txtAdvEvaluation');
uicontrol('parent',hpanelAdvComp,'Units','normalized','style','pushbutton','backgroundcolor',pbBackGroundColor,...
    'position',[0.2 0.4 0.6 0.25],'string','Settings...','FontSize',lFontSize,'enable','off','tag','pbMCRunOptions',...
    'FontName',defaultFontName,'callback','multiRunOptions;');
uicontrol('parent',hpanelAdvComp,'Units','normalized','style','pushbutton','backgroundcolor',pbBackGroundColor,...
    'position',[0.2 0.05 0.6 0.25],'string','Compare','FontSize',lFontSize,'enable','off','tag','pbMCRun',...
    'FontName',defaultFontName,'callback','MCTDRun;');
%% Advanced comparison


%% Quit
hpbExit=uicontrol('parent',hmainTDALAB,'Units','normalized','string','Exit','backgroundcolor',pbBackGroundColor,...
    'style','pushbutton','position',[1-nBorderWidth-0.1 0.02 0.1 0.04],'fontsize',lFontSize,...
    'callback',@TDALABexit);

%% End
refresh(hmainTDALAB);
movegui(hmainTDALAB,'center');
set(hmainTDALAB,'HandleVisibility','callback','visible','on');


%% some callbacks
    function setNonegative(hObject,eventdata)
        tdalabStatus.nonnegativity=(get(hObject,'value')==1);
        updateAlgList;
    end
    function isFullTensor(hObject,eventdata)
        tdalabStatus.fullTensor=(get(hObject,'value')==1);
        if ~tdalabStatus.fullTensor
            cleartemp;
            tdalabStatus.noiseSNR=inf;
            tdalabStatus.noiseType='';
        end
        updateUI;
    end
    function setAlg(hObject,eventdata,handles)
        if ~isempty(tdalabStatus.validAlgs)
            selected=get(hObject,'value');
            tdalabStatus.algIndex=tdalabStatus.validAlgs(selected);
            tdalabStatus.decomposed=false;
            set(hpbOpts,'enable','on');
            set(hgbExtraConstraints,'SelectedObject',[]);
        else
            set(hpbOpts,'enable','off');
        end
    end
    function succeed=selInputFormat(hObject,eventdata)
        ni=get(hObject,'value');
        inputFormat=get(hObject,'string');
        if strcmpi(inputFormat{ni},'ktensor')&&strcmp(tdalabStatus.inputType,'ttensor')
            if all(NumOfComp==NumOfComp(1))
                Y=ktensor(ones(NumOfComp(1),1),Y.U);
                tdalabStatus.inputType='ktensor';
                cleartemp;
                succeed=true;                             
                set(findobj('tag','pmG'),'value',1);
            else
                errordlg('In CP all the mode matrices must have the same column numbers.','Error','modal');
                set(hObject,'value',find(strcmpi(defstr('input'),'ttensor'),1));
                succeed=false;
            end
        elseif strcmpi(inputFormat{ni},'ttensor') && strcmpi(tdalabStatus.inputType,'ktensor')
            subs=(1:NumOfComp(1))';subs=subs(:,ones(1,NumOfMode));
            G=sptensor(subs,Y.lambda);
            Y=ttensor(tensor(G),Y.U);
            tdalabStatus.inputType='ttensor';
            cleartemp;
            succeed=true;    
            set(findobj('tag','pmG'),'value',1);
        end
        updateUI;
    end

    function setNumOfComp(hObject,eventdata)
        if isempty(tdalabStatus.inputType)
            set(hObject,'value',false);
            errordlg('No input. Please load a tensor first.','No input','modal');
            return;
        end
        if get(hObject,'value')
            if ~tdalabStatus.fullTensor
                errordlg('Please use full tensor mode by checking the box [Full tensor].','Error settings','modal');
                set(hObject,'value',false);
                return;
            end
            if any(NumOfComp>1)
                choice = questdlg(strcat('Original number of components is [',...
                    num2str(NumOfComp),']. Continue and update it?'), ...
                        'Please confirm ...','Continue','Cancel','Cancel');
                if strcmp(choice,'Cancel')
                    set(hObject,'value',false);
                    return;
                else
                    pause(0.01);
                end
            end
            tdalab('hide');
            commandwindow;
            fprintf('[TDALAB] Begin to estimate the number of components ...\n');
            if exist(TEMPFILE,'file')~=0
                load(TEMPFILE,'Ynoise','SNR');
            else
                Ynoise=[];SNR=rand;
            end
            if isempty(Ynoise)||(SNR~=tdalabStatus.noiseSNR)
                    if isinf(tdalabStatus.noiseSNR)
                        fprintf('[TDALAB] Generating the full tensor. This may cost a few minutes depending on the problem size. Please wait ...\n');
                        Ynoise=tensor(Y);
                        SNR=tdalabStatus.noiseSNR;
                        save(TEMPFILE,'Ynoise','-v7.3','SNR');
                    else
                        h=errordlg('Noise data not found. Please check your noise settings.','Error settings','modal');
                        uiwait(h);
                        tdalab;
                        return;
                    end
            end
            for n=1:NumOfMode
                [NumOfComp(n) gs{n} es{n}]=moderankdect(Ynoise,n);
            end
            NumOfComp=NOCanalyze(gs,es);
            if strcmpi(tdalabStatus.model,'CP')
                NumOfComp=repmat(max(NumOfComp),1,NumOfMode);
            end
            fprintf('[TDALAB] The estimated number of components (mode-rank): %s.\n',num2str(NumOfComp));
        else
            if strcmpi(tdalabStatus.inputType,'ktensor')
                NumOfComp=numel(Y.lambda);
            elseif strcmpi(tdalabStatus.inputType,'ttensor')
                NumOfComp=size(Y.core);
            else
                NumOfComp=zeros(1,NumOfMode);
            end
        end
        tdalab;
        updateUI;
    end

    function selTDModel(hObject,eventdata)
        oldmodel=tdalabStatus.model;
        tdalabStatus.model=TDModelstr{get(hObject,'value')};
        if strcmpi(tdalabStatus.model,'Tucker')
            %% update input: ktensor-->ttensor
            if strcmpi(tdalabStatus.inputType,'ktensor')
                hinput=findobj('tag','pmInputFormat');
                set(hinput,'value',find(strcmpi(inputFormat,'ttensor'),1));
                selInputFormat(hinput,[]);
            end
        elseif strcmpi(tdalabStatus.model,'CP')&&strcmpi(tdalabStatus.inputType,'ttensor')
                hinput=findobj('tag','pmInputFormat');
                set(hinput,'value',find(strcmpi(inputFormat,'ktensor'),1));   
                succeed=selInputFormat(hinput,[]);
                if ~succeed
                    tdalabStatus.model=oldmodel;
                    set(hObject,'value',find(strcmpi(TDModelstr,tdalabStatus.model),1));
                end
        end
        updateUI;
    end
    function FigureExport(hObject,eventdata,handles)
        exportsetupdlg(hmainTDALAB);
    end
    function saveYcapOnly(hObject,eventdata,handles)
        if ~tdalabStatus.decomposed
            errordlg('No results found. Please decompose a tensor first.','No data found');
            return;
        end
        str=get(hmainTDALAB,'name');
        str=regexp(str,'::','split');
        str=strtrim(strrep(str{1},'.mat',''));
        defaultFileName=horzcat(str,'_',date,'.mat');
        [FileName,PathName]=uiputfile('*.mat',horzcat(TDALABHOME,filesep,'userdata'),defaultFileName);
        if ~isequal(FileName,0)&&~isequal(PathName,0)
            save(horzcat(PathName,FileName),'Ycap','-v7.3');
        end
    end

    function gotoCCW(hObject,eventdata,handles)
    commandwindow;
    disp('[TDALAB] Type or click <a href = "matlab: tdalab">tdalab</a> to return.')
    end
    function closeFigs(hObject,eventdata) 
        delete(findobj('-regexp','tag','\<tdvfig.*'));
        delete(findobj('-regexp','tag','\<mtdvfig.*'));
    end
    function TDALABexit(hObject,eventdata)
        choice = questdlg(strcat('Do you want to quit TDALAB now?'), ...
            'Close TDALAB ...', ...
            'Yes','Save Path & Close','No','No');
        % Handle response
        switch choice
            case 'Yes'
                oldpath=cd(TDALABHOME);
                closeFigs;
                if ~exist('hmainTDALAB')||~ishandle(hmainTDALAB)
                    hall=allchild(0);
                    hmainTDALAB=findall(hall,'tag','TDALAB');
                end
                cleartemp;
                
                rmpath(TDALABPATH);
                
                delete(hmainTDALAB);
                cd(oldpath);
                clear;
                warning('on');
                return;
            case 'Save Path & Close'
                oldpath=cd(TDALABHOME);
                closeFigs;
                if ~exist('hmainTDALAB')||~ishandle(hmainTDALAB)
                    hall=allchild(0);
                    hmainTDALAB=findall(hall,'tag','TDALAB');
                end
                cleartemp;
                
                delete(hmainTDALAB);
                cd(oldpath);
                clear;
                warning('on');
                return;
            case 'No'
                return;
        end
    end
end