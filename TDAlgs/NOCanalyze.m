function NumOfComp=NOCanalyze(gs,es)
global XSpread YSpread backGroundColor;
global defaultFontName sFontSize lFontSize;
global NumOfMode;
NumOfComp=[];

status='eigs';
firstK=3;
ind=cell(1,NumOfMode);
for n=1:NumOfMode
    if ~isempty(gs{n})
        [p ind{n}]=sort(gs{n},'ascend');
        T=min(10,length(gs{n}));
        gs{n}=p(1:T);
        ind{n}=ind{n}(1:T);
        MaxDisp(n)=min(max(ind{n}(1:3))+T,50);
        es{n}=es{n}(1:min(MaxDisp(n),length(es{n})));
    else
        ind{n}=1:length(es{n});
    end
end


%% defin mainWindow
hfigNOC=figure('Units','normalized','Resize','on','toolbar','none','menu','none',...
    'OuterPosition',[0.1 0.1 .75 .5],'CloseRequestFcn',@ExitNOC,'WindowStyle','modal',...
    'tag','NOCanalyze','Name','SORTE','numbertitle','off','color',backGroundColor,'visible','off');
set(hfigNOC,'units','characters');
pos=get(hfigNOC,'Position');

textH=50;

if (pos(3)-NumOfMode-1)/NumOfMode>pos(4)-textH
    pos(4)=(pos(3)-NumOfMode-1)/NumOfMode+textH;
else
    pos(3)=(pos(4)-23)*23+NumOfMode+1;
end


colormap('default');
bkcolor=get(hfigNOC,'Color');

%% Menu
h=8/pos(4);
uicontrol('parent',hfigNOC,'Units','normalized','fontname',defaultFontName,'position',[0.01 1-h,0.98,h*.8],...
    'style','text','string','SORTE: Second  ORder  sTatistic  of  the  Eigenvalues','backgroundcolor',bkcolor,'fontunits','points','fontsize',14);

space=2/pos(4);
base=1-h-space;
h=4/pos(4);
hgbShowValues=uibuttongroup('parent',hfigNOC,'Units','normalized','fontunits','normalized','fontsize',.8,...
    'position',[0.01 base-h 0.48 h],'backgroundcolor',bkcolor,...
    'bordertype','line','borderwidth',0,'visible','on','SelectionChangeFcn',@setShowValues);
hrbEigs=uicontrol('parent',hgbShowValues,'Units','normalized','fontunits','points','fontsize',11,'FontName',defaultFontName,...
    'position',[0.01 0.05 0.48 0.9],'background',bkcolor,'style','radiobutton','tag','rbEigs','string','EigenValues');
hrbGaps=uicontrol('parent',hgbShowValues,'Units','normalized','fontunits','points','fontsize',11,'FontName',defaultFontName,...
    'position',[0.51 0.05 0.48 0.9],'background',bkcolor,'style','radiobutton','tag','rbEigs','string','GAPs');
set(hgbShowValues,'SelectedObject',hrbEigs);

%% plots
colW=(0.99-NumOfMode*0.01)/NumOfMode;
base=base-h-space;
h=(pos(4)-textH)/pos(4);
for n=1:NumOfMode
    labelspace=ceil(length(gs{n})/10);
    hax(n)=axes('parent',hfigNOC,'units','normalized','outerposition',[0.01+(n-1)*colW base-h colW h]);
    plot(es{n},'*-');
    set(gca,'XTick',1:length(es{n}),'XTickLabel',1:labelspace:length(es{n}));
    axis tight;
    xlabel('Index (NumOfComp)');
    ylabel('Eigenvalues');
    grid on;
    title(horzcat('Mode ',num2str(n)));
end

%% mode detect
base=base-h-space;
h=14/pos(4)+2*space;
hk=0.98/firstK;
for n=1:NumOfMode
    hbgM(n)=uibuttongroup('parent',hfigNOC,'Units','normalized','fontsize',lFontSize,'userdata',n,...
        'position',[0.01+(n-1)*colW base-h colW-0.01 h],'backgroundcolor',bkcolor,'visible','on','title',horzcat('Mode ',num2str(n)));
    
    for k=1:firstK
        uicontrol('parent',hbgM(n),'Units','normalized','position',[0.01 1-hk*k 0.98 hk],'fontname',defaultFontName,...
               'fontsize',lFontSize,'string',num2str(ind{n}(k)),'style','radiobutton','backgroundcolor',bkcolor);
    end
end

base=base-h-space;
h=8/pos(4);
hexit=uicontrol('parent',hfigNOC,'Units','normalized','fontname',defaultFontName,'fontsize',lFontSize,...
    'position',[0.4 0.02 0.2 h],'style','pushbutton','string','Confirm','callback',@ExitNOC);
movegui(hfigNOC,'center');
set(hfigNOC,'visible','on');
waitfor(hexit);


    function setShowValues(source,eventdata)
    currentConstr=get(eventdata.NewValue,'String');
    if strcmpi(currentConstr,'GAPs')
        for n=1:NumOfMode
            labelspace=ceil(length(gs{n})/10);
            if ~isempty(gs{n})
                cla(hax(n));
                axes(hax(n));
                plot(gs{n},'*-b');
                set(gca,'XTick',1:labelspace:length(gs{n}));
                ylabel('GAPs');
            end
            axis tight;
            xlabel('Index (NumOfComp)');
            grid on;
            title(horzcat('Mode ',num2str(n)));
        end
    else
        for n=1:NumOfMode
            labelspace=ceil(length(es{n})/10);
            cla(hax(n));
            axes(hax(n));
            plot(es{n},'*-b');
            set(gca,'XTick',1:labelspace:length(es{n}));
            axis tight;
            xlabel('Index (NumOfComp)');
            ylabel('Eigenvalues');
            grid on;
            title(horzcat('Mode ',num2str(n)));
        end
    end
    end

    function setNumOfComp(source,eventdata)   
        NumOfComp(get(source,'userdata'))=str2double(get(eventdata.NewValue,'String'));
    end
        
    function ExitNOC(hObject,eventdata,handle)
        for n=1:NumOfMode
            NumOfComp(n)=str2double(get(get(hbgM(n),'SelectedObject'),'string'));
        end
        delete(hfigNOC);
    end

end