function runclustering
global XSpread YSpread ScreenWidth ScreenHeight defaultCtrlHeight;
global TDALABHOME backGroundColor NumOfMode tdalabStatus;
global Y Ycap;

h=findobj('tag','figTensorClss');
if ~isempty(h)
    figure(h);
    return;
end

if isempty(tdalabStatus.inputType)||isempty(NumOfMode)
    errordlg('Load data first.','No data found');
    return;
end

drm_flag=false;
drm_opts=struct();
kmeans_opts=struct();
gnd_label=[];
est_label=[];
Feas=[];
CPerformance=[];

hmode=[];
hselData=[];
hcore=[];

kmeansalg(1).name='kmeans';kmeansalg(1).details='K-means clustering';
kmeansPList{1}=struct('k',[],'distance','sqEuclidean|cityblock|cosine|correlation|Hamming','emptyaction','error|drop|singleton',...
    'onlinephase','on|off','replicates',5,'start','sample|uniform|cluster','Display','off|final|iter');
kmeansTPList{1}=struct('k',paraType.paraDouble,'distance',paraType.paraMenu,'emptyaction',paraType.paraMenu,...
    'onlinephase',paraType.paraMenu,'replicates',paraType.paraDouble,'start',paraType.paraMenu,'Display',paraType.paraMenu);
kmeansDefault=struct('k',[],'distance','sqEuclidean','emptyaction','error',...
    'onlinephase','on','replicates',5,'start','sample','Display','off');


tdalab('hide');

mainHeight=25*YSpread;mainWidth=105*XSpread;

hcluster=figure('Units','Characters','Resize','on','toolbar','none','menu','none','color',backGroundColor,...
    'OuterPosition',[ScreenWidth/2-mainWidth/2,ScreenHeight-mainHeight-defaultCtrlHeight,mainWidth,mainHeight],...
    'tag','figTensorClss','Name','Tensor Clustering Analysis','numbertitle','off','visible','off','CloseRequestFcn',@clusterExit);
colormap('default');
bkcolor=get(hcluster,'Color');

upbase=0.95;
downspace=0.01;
ctrlHS=(upbase-downspace)/6.5;
ctrlH=ctrlHS*.7;
textH=0.5;
left=0.1;

%% title
uicontrol('parent',hcluster,'Units','normalized','position',[.1 upbase-1.25*ctrlHS .8 1.5*ctrlH],'string','Tensor Clustering Analysis','style','text','fontunits','normalized','FontSize',.5,...
    'BackgroundColor',bkcolor,'fontweight','bold');

%% data
modestr=num2cell(1:NumOfMode);
upbase=upbase-2.5*ctrlHS;
uicontrol('parent',hcluster,'Units','normalized','position',[left upbase .12 ctrlH],'string','Mode:','style','text','fontunits','normalized','FontSize',textH,...
    'BackgroundColor',bkcolor,'horizontalalign','left');
hmode=uicontrol('parent',hcluster,'Units','normalized','position',[left+.125 upbase+ctrlH*.5 .2 ctrlH*.8],'string',modestr,'style','popupmenu','fontunits','normalized','FontSize',textH,'value',NumOfMode);
hselData=uicontrol('parent',hcluster,'Units','normalized','position',[left+.4 upbase .3 ctrlH],'string','Decomposed data','style','checkbox','fontunits','normalized','FontSize',textH,...
    'backgroundcolor',bkcolor,'max',true,'min',false,'value',true);
hcore=uicontrol('parent',hcluster,'Units','normalized','position',[left+.71 upbase .25 ctrlH],'string','with core','style','checkbox','fontunits','normalized','FontSize',textH,...
    'backgroundcolor',bkcolor,'max',true,'min',false,'value',true);

if tdalabStatus.decomposed
    set([hselData hcore],'enable','on');
else %% using the original data
    set(hselData,'value',true,'enable','off');
    if strcmpi(class(Y),'tensor')
        set(hcore,'enable','off');
    else
        set(hcore,'enable','on');
    end
end
        

%% dimensionality reduction
upbase=upbase-ctrlHS;
hdrm=uicontrol('parent',hcluster,'Units','normalized','position',[left upbase .4 ctrlH],'string','Dimensionality reduction','style','checkbox','fontunits','normalized','FontSize',textH,...
    'backgroundcolor',bkcolor,'callback',@setdimrm);

%% set kmeans
upbase=upbase-1.25*ctrlHS;
uicontrol('parent',hcluster,'Units','normalized','position',[left upbase .4 ctrlH],'string','KMeans configuration ...','style','pushbutton','fontunits','normalized','FontSize',textH,...
    'callback',@setkmeans);

%% load true lable
uicontrol('parent',hcluster,'Units','normalized','position',[.55 upbase .4 ctrlH],'string','Load true label ... ','style','pushbutton','fontunits','normalized','FontSize',textH,...
    'callback',@loadTrueLabel);

%% run -vis -save -exit
upbase=upbase-1.5*ctrlHS;
ctrlw=0.2;
uicontrol('parent',hcluster,'Units','normalized','position',[.1 upbase ctrlw ctrlH],'string','Run','style','pushbutton','fontunits','normalized','FontSize',textH,...
    'callback',@runKMeans);
uicontrol('parent',hcluster,'Units','normalized','position',[.4 upbase ctrlw ctrlH],'string','Save','style','pushbutton','fontunits','normalized','FontSize',textH,...
    'callback',@clusterSave);
uicontrol('parent',hcluster,'Units','normalized','position',[.7 upbase ctrlw ctrlH],'string','Exit','style','pushbutton','fontunits','normalized','FontSize',textH,...
    'callback',@clusterExit);


set(hcluster,'visible','on');
movegui(hcluster,'center');

    function runKMeans(hObject,eventdata,handles)
        
        delete(findobj('-regexp','tag','\<tdvfig.*'));
        
        %% get features
        if get(hselData,'value')
            switch tdalabStatus.model
                case 'CP'
                    Feas=Ycap.U{get(hmode,'value')};
                    if get(hcore,'value')
                        Feas=bsxfun(@times,Feas,Ycap.lambda');
                    end
                case 'Tucker'
                    Feas=Ycap.U{get(hmode,'value')};
                    if get(hcore,'value')
                        Feas=ttm(Ycap.core,Feas,get(hmode,'value'));
                        Feas=tenmat(Feas,get(hmode,'value'));
                        Feas=double(Feas);
                    end
                otherwise
                    errordlg('Unsupported tensor model.');
                    return;
            end
        else % not decomposed
            switch class(Y)
                case 'CP'
                    Feas=Y.U{get(hmode,'value')};
                    if get(hcore,'value')
                        Feas=bsxfun(@times,Feas,Y.lambda);
                    end
                case 'Tucker'
                    Feas=Y.U{get(hmode,'value')};
                    if get(hcore,'value')
                        Feas=ttm(Y.core,Feas,get(hmode,'value'));
                        Feas=tenmat(Feas,get(hmode,'value'));
                        Feas=double(Feas);
                    end
                otherwise
                    fprintf(2,'The tensor has not been decomposed yet. Original data will be used.\n');
                    Feas=double(tenmat(Y,get(hmode,'value')));
            end           
        end
        
        set(hcluster,'visible','off');
        commandwindow;
        
        %% k-means parameters
        [k,distance,emptyaction,onlinephase,replicates,kmstart,Display]=scanparam(kmeansDefault,kmeans_opts);
        if isempty(k)
            if isempty(gnd_label)
                h=errordlg('Number of clusters must be specified.','Key parameters missing');
                uiwait(h);
                set(hcluster,'visible','on');
                return;
            else
                k=numel(unique(gnd_label));
            end
        end     
        %% dim reduction
        if drm_flag
            disp('Dimensionality reduction is running ... ');
            disp(drm_opts)
            try
                Feas=PMFcompute_mapping(Feas,drm_opts);
            catch ME
                for id=numel(ME.stack):-1:1
                    disp(ME.stack(id));
                end
            end                
        end
        

        disp('K-Means is running. Please wait ...');
        kmopts= statset('Display',Display);   
        try
            est_label=kmeans(Feas,k,'distance',distance,'emptyaction',emptyaction,'onlinephase',onlinephase,'replicates',replicates,'start',kmstart,'options',kmopts);
        catch ME
            for id=numel(ME.stack):-1:1
                disp(ME.stack(id));
            end
        end
        
        if ~isempty(gnd_label)
            commandwindow;
           
            disp('Calculating the clustering accuracy ... ');
            fprintf('Clustering Accuracy: %3.1f%%\n',accuracy(gnd_label,est_label));
        end
        
        if size(Feas,2)<=3
            hcv=visualizeClusters(Feas,est_label,gnd_label);
            set([hcv hcluster],'units','normalized');
            posL=get(hcluster,'outerposition');
            posR=get(hcv,'outerposition');
            rs=posL(4)/posR(4);
            posR(3:4)=posR(3:4)*rs;
            left=(1-posL(3)-posR(3))/2;
            set(hcluster,'outerposition',[left posL(2:4)]);
            set(hcv,'outerposition',[left+posL(3) posL(2) posR(3:4)]);
        end
        
        set(hcluster,'visible','on');
    end

    function clusterSave(hObject,eventdata,handles)
        [FileName FilePath]=uiputfile('*.mat',horzcat(TDALABHOME,filesep,'userdata',filesep));
        if  ~isequal(FileName,0) && ~isequal(FilePath,0)
            if isempty(CPerformance)
                save(fullfile(FilePath,FileName),'est_label','Feas');
            else
                save(fullfile(FilePath,FileName),'est_label','Feas','CPerformance');
            end
        end
    end

    function loadTrueLabel(hObject,eventdata,handles)
        [FileName FilePath]=uigetfile('*.mat','Select a file including the true label');
        if ~isequal(FileName,0)
            ls=load(fullfile(FilePath,FileName),'label');
            if isfield(ls,'label')
                gnd_label=reshape(ls.label,[],1);
                kmeansDefault.k=numel(unique(gnd_label));
            end                
        end
        if ~isempty(gnd_label)
            set(hObject,'string','True label available.');
        else
            set(hObject,'string','Load true label ...');
        end
    end

    function setkmeans(hObject,eventdata,handles)
        kmeans_opts=guiSetOpts(kmeansalg,kmeansTPList,kmeansPList,kmeansDefault);
        if ~isempty(kmeans_opts)
            kmeans_opts=kmeans_opts{1};
        else
            kmeans_opts=struct();
        end
    end

    function setdimrm(hObject,eventdata,handles)
        if get(hObject,'value')
            drm_flag=true;
            %% set parameters
            alg(1).name='PMFcompute_mapping';alg(1).details='Dimensionality reduction';
            algPList{1}=struct('method','tSNE|PCA|ISOMAP|LLE','no_dims',3,'NNeighbors',12,'init_dims',30);
            algTPList{1}=struct('method',paraType.paraMenu,'no_dims',paraType.paraDouble,'NNeighbors',paraType.paraDouble,'init_dims',paraType.paraDouble);
            drm_default_opts=struct('method','tSNE','no_dims',2);
            drm_opts=guiSetOpts(alg,algTPList,algPList);
            if ~isempty(drm_opts)
                drm_opts=drm_opts{1};
            else
                drm_opts=drm_default_opts;
            end
        else
            drm_flag=false;
            drm_opts=struct();
        end
    end
    function clusterExit(hObject,eventdata,handles)
        delete(findobj('-regexp','tag','\<tdvfig.*'));
        delete(hcluster);
        tdalab('show');
    end
end