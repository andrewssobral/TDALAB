function callTuckerDA
tdalab('hide');
global XSpread YSpread ScreenWidth ScreenHeight defaultCtrlHeight;
global sFontSize lFontSize;
global TDALABHOME backGroundColor;

h=findobj('tag','figtensorDA');
if h
    figure(h);
    return;
end

data=[];
hpara_NOC=[];
hpara_maxit=[];
hpara_tol=[];
hpara_eig=[];
hpara_classifier=[];
hpara_classify_rule=[];
hpara_classify_dist=[];

evaluationMode=false;
Accuracy=-1;
label_hat=[];
CPerformance=[];
sample_fea=[];
training_fea=[];
tensorDA_paras=struct();
hdatainfo=[];

datafile=[];

lda_rules={'linear','diaglinear','quadratic','diagquadratic','mahalanobis'};
knn_rules={'nearest','random','consensus'};
distance={'euclidean','cityblock','cosine','correlation','hamming'};

mainHeight=45*YSpread;mainWidth=110*XSpread;
ctrlHeight=.07;

htensorDA=figure('Units','Characters','Resize','on','toolbar','none','menu','none','color',backGroundColor,...
    'OuterPosition',[ScreenWidth/2-mainWidth/2,ScreenHeight-mainHeight-defaultCtrlHeight,mainWidth,mainHeight],...
    'tag','figtensorDA','Name','Tensor Discriminant Analysis','numbertitle','off','visible','off','CloseRequestFcn',@tensorDAexit);
colormap('default');
bkcolor=get(htensorDA,'Color');
uicontrol('parent',htensorDA,'Units','normalized','position',[.1 .85 .8 .1],'string','Tucker Discriminant Analysis','style','text','FontSize',lFontSize*1.5,...
    'BackgroundColor',bkcolor,'fontweight','bold');

uicontrol('parent',htensorDA,'Units','normalized','position',[.01 .81-ctrlHeight 0.36 ctrlHeight],'string','Load data','FontSize',lFontSize,...
    'style','pushbutton','CallBack',@loadTensorDAData,'BackgroundColor',bkcolor);
hdatainfo=uicontrol('parent',htensorDA,'Units','normalized','position',[.01 .125 0.36 0.58],'string','-- Load data first --','FontSize',lFontSize,...
    'style','edit','min',1,'max',4,'background',bkcolor,'enable','inactive','horizontalalign','left',...
    'fontunits','point','fontsize',11);

hpPara=uipanel('parent',htensorDA,'Units','normalized','position',[.38 .125 0.6 .7],'title','Parameters','FontSize',lFontSize,...
    'FontUnits','points','fontsize',sFontSize,'foregroundcolor','k','backgroundcolor',bkcolor,'borderwidth',1,'highlightcolor','k',...
    'backgroundcolor',bkcolor,'BorderType','line');

uicontrol('parent',htensorDA,'Units','normalized','position',[.1 .01 .2 ctrlHeight],'string','Run', 'FontSize',lFontSize,...
    'style','pushbutton','callback',@runTensorDA);
uicontrol('parent',htensorDA,'Units','normalized','position',[.4 .01 .2 ctrlHeight],'string','Save results', 'FontSize',lFontSize,...
    'style','pushbutton','callback',@saveTensorDA);
uicontrol('parent',htensorDA,'Units','normalized','position',[.7 .01 .2 ctrlHeight],'string','Quit', 'FontSize',lFontSize,...
    'style','pushbutton','callback',@tensorDAexit);

%% parameters
ctrlHeight=.1;
sctrlHeight=.3/8+ctrlHeight;

% label
uicontrol('parent',hpPara,'Units','normalized','position',[.1 1-sctrlHeight 0.38 ctrlHeight],'string','Dim. of Features:','FontSize',lFontSize,...
    'backgroundcolor',bkcolor,'horizontalalignment','right','style','text');
uicontrol('parent',hpPara,'Units','normalized','position',[.1 1-2*sctrlHeight 0.38 ctrlHeight],'string','Max. iteration:','FontSize',lFontSize,...
    'backgroundcolor',bkcolor,'horizontalalignment','right','style','text');
uicontrol('parent',hpPara,'Units','normalized','position',[.1 1-3*sctrlHeight 0.38 ctrlHeight],'string','tol:','FontSize',lFontSize,...
    'backgroundcolor',bkcolor,'horizontalalignment','right','style','text');
uicontrol('parent',hpPara,'Units','normalized','position',[.1 1-4*sctrlHeight 0.38 ctrlHeight],'string','Eig Mode:','FontSize',lFontSize,...
    'backgroundcolor',bkcolor,'horizontalalignment','right','style','text');
uicontrol('parent',hpPara,'Units','normalized','position',[.1 1-5*sctrlHeight 0.38 ctrlHeight],'string','Classifier:','FontSize',lFontSize,...
    'backgroundcolor',bkcolor,'horizontalalignment','right','style','text');
uicontrol('parent',hpPara,'Units','normalized','position',[.1 1-6*sctrlHeight 0.38 ctrlHeight],'string','Classify rule:','FontSize',lFontSize,...
    'backgroundcolor',bkcolor,'horizontalalignment','right','style','text');
uicontrol('parent',hpPara,'Units','normalized','position',[.1 1-7*sctrlHeight 0.38 ctrlHeight],'string','Distance:','FontSize',lFontSize,...
    'backgroundcolor',bkcolor,'horizontalalignment','right','style','text');
% value
% ctrlHeight=.08;
hpara_NOC=uicontrol('parent',hpPara,'Units','normalized','position',[.5 1.025-sctrlHeight 0.45 ctrlHeight],'string','Projected dim. of each feature mode','FontSize',lFontSize,...
    'backgroundcolor',[1 1 1],'horizontalalignment','left','style','edit','ButtonDownFcn',@clearField,'callback',@setNumOfComp);
hpara_maxit=uicontrol('parent',hpPara,'Units','normalized','position',[.5 1.025-2*sctrlHeight 0.45 ctrlHeight],'string','200','FontSize',lFontSize,...
    'backgroundcolor',[1 1 1],'horizontalalignment','left','style','edit','ButtonDownFcn',@clearField);
hpara_tol=uicontrol('parent',hpPara,'Units','normalized','position',[.5 1.025-3*sctrlHeight 0.45 ctrlHeight],'string','0.001','FontSize',lFontSize,...
    'backgroundcolor',[1 1 1],'horizontalalignment','left','style','edit','ButtonDownFcn',@clearField);
hpara_eig=uicontrol('parent',hpPara,'Units','normalized','position',[.5 1.025-4*sctrlHeight 0.45 ctrlHeight],'string',{'ratio','diff'},'FontSize',lFontSize,...
    'backgroundcolor',[1 1 1],'horizontalalignment','left','style','popupmenu');
hpara_classifier=uicontrol('parent',hpPara,'Units','normalized','position',[.5 1.025-5*sctrlHeight 0.45 ctrlHeight],'string',{'LDA','kNN'},'FontSize',lFontSize,...
    'backgroundcolor',[1 1 1],'horizontalalignment','left','style','popupmenu','callBack',@setClassifier);
hpara_classify_rule=uicontrol('parent',hpPara,'Units','normalized','position',[.5 1.025-6*sctrlHeight 0.45 ctrlHeight],'string',lda_rules,'FontSize',lFontSize,...
    'backgroundcolor',[1 1 1],'horizontalalignment','left','style','popupmenu');
hpara_classify_dist=uicontrol('parent',hpPara,'Units','normalized','position',[.5 1.025-7*sctrlHeight 0.45 ctrlHeight],'string',distance,'FontSize',lFontSize,...
    'backgroundcolor',[1 1 1],'horizontalalignment','left','style','popupmenu','enable','off');

movegui(htensorDA,'center');
set(htensorDA,'visible','on');

    function runTensorDA(hObject,eventdata,handles)
        if isempty(data)
            errordlg('Load data first.','Error','modal');
            return;
        end
        strs=get(hpara_classifier,'string');
        classifier=strs{get(hpara_classifier,'value')};
        classify_distance=distance{get(hpara_classify_dist,'value')};
        if strcmpi(classifier,'lda')
            classify_rule=lda_rules{get(hpara_classify_rule,'value')};
        elseif strcmpi(classifier,'knn')
            if exist('knnclassify','file')
                classify_rule=knn_rules{get(hpara_classify_rule,'value')};
            else
                errordlg('knnclassify does not exist.','File not found','modal');
                return;
            end
        end
        NumOfComp=str2num(get(hpara_NOC,'string'));
        maxit=str2double(get(hpara_maxit,'string'));
        tol=str2double(get(hpara_tol,'string'));
        strs=get(hpara_eig,'string');
        eigmode=strs{get(hpara_eig,'value')};
        tensorDA_paras=struct('NumOfComp',NumOfComp,'maxiter',maxit,'tol',tol,'Eigs',eigmode,'classifier',classifier,'distance',classify_distance,'classify_rule',classify_rule);
        commandwindow;
        disp('Tensor Discriminant Analysis is running ...');
        [label_hat sample_fea training_fea]=tuckerDA(data.sample,data.training,data.group,tensorDA_paras);
        disp('Completed.');
        if evaluationMode
            if exist('classperf','file')
                CPerformance=classperf(data.label,label_hat);
                disp('Detailed classification performance:');
                CPerformance
            end
            fprintf('\n==================================\n');
            Accuracy=sum(label_hat==data.label)./numel(data.label);
            fprintf('Accuracy: %3.1f%%\n',100*Accuracy);
            fprintf('==================================\n');
        end
        figure(htensorDA);
    end

    function setNumOfComp(hObject,eventdata,handles)
        if ~isempty(data)
            datasize=size(data.training);
            noc=round(str2num(get(hpara_NOC,'string')));
            NOMode=numel(datasize);
            if numel(noc)~=NOMode-1
                noc(end+1:NOMode-1)=noc(end);
                noc=noc(1:NOMode-1);
            end
            set(hObject,'string',num2str(noc));
        end
    end

    function saveTensorDA(hObject,eventdata,handles)
        if isempty(label_hat)
            errordlg('Run Tucker Discriminant Analysis first.','No results available','modal');
            return;
        end
        [FileName PathName]=uiputfile('*.mat','Save results as ...',horzcat(TDALABHOME,filesep,'userdata',filesep,'TDA-',datafile,'-',date,'.mat'));
        if ~isequal(FileName,0)&&~isequal(PathName,0)
            varlist={'label_hat','tensorDA_paras','datafile','Accuracy','sample_fea','training_fea'};
            if ~isempty(CPerformance)
                varlist{end+1}='CPerformance';
            end
            timestamp=horzcat('Saved at ',datestr(now));
            varlist{end+1}='timestamp';
            
            %% construct values
            for n=1:numel(varlist)
                if ~isempty(varlist{n})
                    TDA.(varlist{n})=eval(varlist{n});
                end
            end
            
            save(horzcat(PathName,FileName),'TDA');
        end % to save
    end

    function setClassifier(hObject,eventdata,handles)
        strs=get(hObject,'string');
        if strcmpi(strs{get(hObject,'value')},'lda')
            set(hpara_classify_rule,'string',lda_rules,'value',1);
            set(hpara_classify_dist,'enable','off');
        elseif strcmpi(strs{get(hObject,'value')},'knn')
            set(hpara_classify_rule,'string',knn_rules,'value',1);
            set(hpara_classify_dist,'enable','on');
        end
    end
    function clearField(hObject,eventdata,handles)
            set(hObject,'string','');
    end

    function loadTensorDAData(hObject,eventdata,handles)
        [datafile,PathName]=uigetfile('*.mat','Load data for discriminant analysis ...',horzcat(TDALABHOME,filesep,'benchmark',filesep,'*.mat'));
        if ~isequal(datafile,0)
            data=load(horzcat(PathName,datafile));
            if isfield(data,'label')
                evaluationMode=true;
            else
                evaluationMode=false;
            end
            if any(isfield(data,{'sample','training','group'})==0)
                errordlg('A valid data file for discriminant analysis must contain the variables ''sample'',''training'',and ''group''. [The variable ''label'' for true label information is optional.]','Error data file.');
                return;
            end
            szS=size(data.sample);
            szT=size(data.training);
            str{1}=horzcat('Dim. of samples: [',num2str(szS(1:end-1)),']');
            str{2}=horzcat('Number of samples: ',num2str(szS(end)),num2str(100*szS(end)/(szS(end)+szT(end)),'  [%3.1f%%]'));
            str{3}=horzcat('Number of training: ',num2str(szT(end)),num2str(100*szT(end)/(szS(end)+szT(end)),' [%3.1f%%]'));
            if isfield(data,'label')
                str{4}='True lables are known.';
            else
                str{4}='True lables are unknown.';
            end
            str{5}='';
            set(hdatainfo,'string',str);
        end
    end

    function tensorDAexit(hObject,eventdata)
        delete(htensorDA);
        tdalab('show');
    end
end