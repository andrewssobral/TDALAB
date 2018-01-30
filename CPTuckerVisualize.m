function CPTuckerVisualize
global tdalabStatus Y Ycap SIRs NumOfMode;
global XSpread YSpread backGroundColor defaultFontName lFontSize;
modePlotstr={'plot3','plot','waterfall','mplot','stem3','ribbon'};
imgPlotstr={'imagesc','contourf','surf'};

if isempty(Y)&&isempty(Ycap)
    errordlg('No tensors for visualization.','No input','modal');
    return;
end


h=findobj('tag','mtdvfigMenuCPTucker');
if ~isempty(h)
    figure(h);
    h=findobj('-regexp','tag','\<tdvfig.*');
    set(h,'visible','on');
    return;
end
delete(findobj('-regexp','tag','\<tdvfig.*'));


hfigVMenu=[];htgWave=[];vstype=[];
currentMode=1;currentPlotFunc='plot3';

hfpmMode=[];hpmSlice=[];
heditXdims=[];heditYdims=[];heditZdims=[];
hpb3DCore=[];drawSlice=false;
%% prepare for plot
ROWS=1;
allPlot=struct('signal',[],'title',['Source signals - mode ' num2str(currentMode)],'tag','tdvfigsap1',...
    'opos',[]);
detailsPlot=struct( 'signal', [], ...
    'standard',0,...
    'x_limit',[],...
    'x_title',[],...
    'num', 1, ...
    'stitle', 'No input', ...
    'sname', 'a', ...  % changed from 'x'
    'title_first_plot', [], ...
    'sposition', [], ...  %[ 0.01 0.05 0.48 0.37 ]
    'sunits', 'normalized', ...
    'tag','tdvfigsdp1',...
    'window_number', 0 );

%% prepare for 3D tensor
Xdims=1;Ydims=2;Zdims=3:NumOfMode;
if strcmpi(tdalabStatus.model,'Tucker')
    corestr='3D Core tensor';
    core3Dplot=struct('signal',[],'name',corestr,'opos',[],'tag','tdvfigsG3D');
    currentImgFunc='imagesc';
    coreslice=struct('title',corestr,'imagefunc',currentImgFunc,'colormap',...
        'Jet','width',[],'cols',[],'stitle','slice','tag','tdvfigsGslice');
    orderChanged=true;
else
    corestr='Rank-1 Terms';
end


width=30*XSpread;
ctrlHeight=2.5*YSpread;
height=5*ctrlHeight+3*YSpread;
hfigVMenu=figure('Units','characters','position',[0 0 width height],'color',backGroundColor,'toolbar','none',...
    'menu','none','Name','Visualization','numbertitle','off','visible','on','tag','mtdvfigMenuCPTucker','windowstyle','normal');
movegui(hfigVMenu,'northwest');


uicontrol('parent',hfigVMenu,'Units','normalized','position',[0.01 0.75 0.98 0.2],'style','text',...
    'string','Visualization','fontname',defaultFontName,'fontsize',lFontSize+2,'backgroundcolor',backGroundColor,...
    'fontweight','bold');

hbgModeExp=uibuttongroup('parent',hfigVMenu,'Units','normalized','position',[0.01 0.18 0.98 0.64],'borderwidth',0,...
    'title','','fontname',defaultFontName,'fontsize',lFontSize,'backgroundcolor',backGroundColor,...
    'highlightcolor','k','BorderType','line','titleposition','centertop','SelectionChangeFcn',@tdvselChanged);
htgWave=uicontrol('parent',hbgModeExp,'Units','normalized','position',[0.01 0.62 0.98 0.3],'style','togglebutton',...
    'string','Waveform','fontname',defaultFontName,'fontsize',lFontSize,'tag','tdvtgtWave');

hsir=uicontrol('parent',hbgModeExp,'Units','normalized','position',[0.01 0.32 0.98 0.3],'style','togglebutton',...
    'string','SIR plot','fontname',defaultFontName,'fontsize',lFontSize,'tag','tdvtgSpec');
if isempty(SIRs)
    set(hsir,'enable','off');
else
    set(hsir,'enable','on');
end

uicontrol('parent',hbgModeExp,'Units','normalized','position',[0.01 0.02 0.98 0.3],'style','togglebutton',...
    'string',corestr,'fontname',defaultFontName,'fontsize',lFontSize,'tag','tdvtglink');

uicontrol('parent',hfigVMenu,'Units','normalized','position',[0.01 0.02 0.98 0.16],'style','pushbutton',...
    'string','TDALAB','fontname',defaultFontName,'fontsize',lFontSize,...
    'callback','tdalab;delete(findobj(''-regexp'',''tag'',''\<tdvfig.*'')); delete(findobj(''-regexp'',''tag'',''\<mtdvfig.*''));');

set(hfigVMenu,'units','characters');
pos=get(hfigVMenu,'outerposition');
smheight=5*ctrlHeight;
hfpmMode=figure('units','characters','outerposition',[0 pos(2)-smheight-0.01 pos(3) smheight],'menu','none','toolbar','none',...
    'numbertitle','off','tag','mtdvfigsmmenu','color',backGroundColor,'visible','on');
uicontrol('parent',hfpmMode,'units','normalized','position',[0.02 0.75 0.96 0.2],'style','text','backgroundcolor',backGroundColor,...
    'Fontname',defaultFontName,'fontsize',lFontSize,'string','Mode matrix');
uicontrol('parent',hfpmMode,'units','normalized','position',[0.02 0.55 0.96 0.2],'style','popupmenu','backgroundcolor','white',...
    'Fontname',defaultFontName,'fontsize',lFontSize,'string',1:NumOfMode,'callback',@setModePlot,'value',1,'tag','tdvpmMode');
uicontrol('parent',hfpmMode,'units','normalized','position',[0.02 0.25 0.96 0.2],'style','text','backgroundcolor',backGroundColor,...
    'Fontname',defaultFontName,'fontsize',lFontSize,'string','Plot function');
uicontrol('parent',hfpmMode,'units','normalized','position',[0.02 0.02 0.96 0.2],'style','popupmenu','backgroundcolor','white',...
    'Fontname',defaultFontName,'fontsize',lFontSize,'string',modePlotstr,'callback',@setMPlotFunc,'value',1,'tag','tdvpmFunc');

                plotUnits='inch';
                oldUnits=get(0,'units');
                set(0,'units',plotUnits);
                WinSize=get(0,'screensize');
                set(0,'units',oldUnits);
                WinSize=WinSize(1,:);
                set(hfigVMenu,'units',plotUnits);
                mpos=get(hfigVMenu,'outerposition');
                wspace=0.01;hspace=0.01;
                top=mpos(2)+mpos(4);
                left=mpos(1)+mpos(3)+wspace;
                swidth=min((WinSize(3)-wspace-mpos(3))/2,4);
                sheight=min((WinSize(4)-hspace*2)/2,4);
                plotPos1=[left top-sheight swidth sheight];
                plotPos2=[left+swidth top-sheight max(swidth,WinSize(3)-left-swidth) sheight];
                plotPos3=[left top-sheight*2 swidth sheight];
                plotPos4=[left+swidth top-sheight*2 max(swidth,WinSize(3)-left-swidth) sheight];


CPTuckerDraw('new');

smheight=8*ctrlHeight;
hfpm3DCore=figure('units','characters','outerposition',[0 pos(2)-smheight-0.01 pos(3) smheight],'menu','none','toolbar','none',...
    'numbertitle','off','tag','mtdvfigsm3DCore','color',backGroundColor,'visible','off','name','3D tensor');
uicontrol('parent',hfpm3DCore,'units','normalized','position',[0.02 0.8 0.96 0.1],'style','text','backgroundcolor',backGroundColor,...
    'Fontname',defaultFontName,'fontsize',lFontSize,'string',horzcat('Tensor order: ',num2str(NumOfMode)),'horizontalalignment','left');
uicontrol('parent',hfpm3DCore,'units','normalized','position',[0.02 0.65 0.45 0.1],'style','text','backgroundcolor',backGroundColor,...
    'Fontname',defaultFontName,'fontsize',lFontSize,'string','Dims in X: ','horizontalalignment','left');
heditXdims=uicontrol('parent',hfpm3DCore,'units','normalized','position',[0.47 0.65 0.51 0.11],'style','edit',...
    'Fontname',defaultFontName,'fontsize',lFontSize,'string',num2str(Xdims),'horizontalalignment','left','callback',@setXdim);
uicontrol('parent',hfpm3DCore,'units','normalized','position',[0.02 0.5 0.45 0.1],'style','text','backgroundcolor',backGroundColor,...
    'Fontname',defaultFontName,'fontsize',lFontSize,'string','Dims in Y: ','horizontalalignment','left');
heditYdims=uicontrol('parent',hfpm3DCore,'units','normalized','position',[0.47 0.5 0.51 0.11],'style','edit',...
    'Fontname',defaultFontName,'fontsize',lFontSize,'string',num2str(Ydims),'horizontalalignment','left','callback',@setYdim);
uicontrol('parent',hfpm3DCore,'units','normalized','position',[0.02 0.35 0.45 0.1],'style','text','backgroundcolor',backGroundColor,...
    'Fontname',defaultFontName,'fontsize',lFontSize,'string','Dims in Z: ','horizontalalignment','left');
heditZdims=uicontrol('parent',hfpm3DCore,'units','normalized','position',[0.47 0.35 0.51 0.11],'style','edit',...
    'Fontname',defaultFontName,'fontsize',lFontSize,'string',num2str(Zdims),'horizontalalignment','left','enable','inactive');
uicontrol('parent',hfpm3DCore,'units','normalized','position',[0.02 0.2 0.35 0.12],'style','checkbox','backgroundcolor',backGroundColor,...
    'Fontname',defaultFontName,'fontsize',lFontSize,'string','slice','horizontalalignment','left','value',false,'max',true,'min',false,'callback',@setSlice);
hpmSlice=uicontrol('parent',hfpm3DCore,'units','normalized','position',[0.38 0.2 0.6 0.12],'style','popupmenu','callback',@setImgFunc,...
    'Fontname',defaultFontName,'fontsize',lFontSize,'string',imgPlotstr,'horizontalalignment','left','value',1,'enable','off');
hpb3DCore=uicontrol('parent',hfpm3DCore,'units','normalized','position',[0.25 0.02 0.5 0.14],'style','pushbutton',...
    'Fontname',defaultFontName,'fontsize',lFontSize,'string','Draw now','callback',@draw3DNow,'enable','on');

    function tdvselChanged(hObject,eventdata)
        if strcmpi(eventdata.EventName,'SelectionChanged')
            delete(findobj('-regexp','tag','\<tdvfigs.*'));
            vstype=get(eventdata.NewValue,'string');
            switch vstype
                case 'Waveform'
                    set(hfpm3DCore,'visible','off');
                    set(hfpmMode,'visible','on');
                    CPTuckerDraw('mode');
                case 'SIR plot'
                    set(hfpm3DCore,'visible','off');
                    set(hfpmMode,'visible','off');
                    showPIs;
                case corestr
                    set(hfpmMode,'visible','off');
                    if strcmpi(tdalabStatus.model,'Tucker')
                        set(hfpm3DCore,'visible','on');
                        delete(findobj('-regexp','tag','\<tdvfig.*'));
                        orderChanged=true;
                        show3Dtensor;
                    else
                        showLambda;
                    end
            end
        end
    end
    function setModePlot(hObject,eventdata)
        currentMode=get(hObject,'value');
        CPTuckerDraw('mode');
    end
    function setMPlotFunc(hObject,eventdata)
        currentPlotFunc=modePlotstr{get(hObject,'value')};
        CPTuckerDraw('func');
    end
    function setImgFunc(hObject,eventdata)
        currentImgFunc=imgPlotstr{get(hObject,'value')};
    end
    function setSlice(hObject,eventdata)
        if get(hObject,'value')
            set(hpmSlice,'enable','on');
            drawSlice=true;
        else
            set(hpmSlice,'enable','off');
            drawSlice=false;
        end
    end
    function setXdim(hObject,eventdata)
        Xdims=intersect(str2num(get(hObject,'string')),1:NumOfMode);
        Ydims=setdiff(str2num(get(heditYdims,'string')),Xdims);
        if isempty(Ydims)
            Ydims=setdiff(1:NumOfMode,Xdims);
            Ydims=Ydims(1);
        end
        set(heditYdims,'string',num2str(Ydims))
        Zdims=setdiff(1:NumOfMode,[Xdims Ydims]);
        set(heditZdims,'string',num2str(Zdims));
        if isempty(Xdims)||isempty(Ydims)||isempty(Zdims)
            set(hpb3DCore,'enable','off');
        else
            set(hpb3DCore,'enable','on');
            orderChanged=true;
        end
    end
    function setYdim(hObject,eventdata)
        Ydims=intersect(str2num(get(hObject,'string')),1:NumOfMode);
        Xdims=setdiff(str2num(get(heditXdims,'string')),Ydims);
        if isempty(Xdims)
            Xdims=setdiff(1:NumOfMode,Ydims);
            Xdims=Xdims(1);
        end
        set(heditYdims,'string',num2str(Ydims))
        Zdims=setdiff(1:NumOfMode,[Xdims Ydims]);
        set(heditZdims,'string',num2str(Zdims));
        if isempty(Xdims)||isempty(Ydims)||isempty(Zdims)
            set(hpb3DCore,'enable','off');
        else
            set(hpb3DCore,'enable','on');
            orderChanged=true;
        end
    end
    function draw3DNow(hObject,eventdata)
        show3Dtensor;
    end

%% other funs
    function CPTuckerDraw(type)
        switch type
            case {'new','mode'}
                ROWS=1;
                if ~isempty(Ycap)&&all(strcmpi(class(Y),{'ktensor','ttensor'})==0)
                    detailsPlot.signal=(Ycap.U{currentMode})';
                    detailsPlot.num=size((Ycap.U{currentMode}),2);
                    detailsPlot.stitle=['Estimated signals - mode ' num2str(currentMode)];
                    detailsPlot.sposition=plotPos1;

                    allPlot.signal=Ycap.U{currentMode}';
                    allPlot.title=['Estimated signals - mode ' num2str(currentMode)];
                    allPlot.opos=plotPos2;
                    
                    allPlot.zlabel=horzcat('{\bf{a}}^{(',num2str(currentMode),')}_{\itj}');
                    allPlot.legend='';
                    for idx=1:size(Ycap.U{currentMode},2)
                        allPlot.legend=horzcat(allPlot.legend,'|{\bf{a}}^{(',num2str(currentMode),')}_',num2str(idx));
                    end                    
                    
                elseif any(strcmpi(class(Y),{'ktensor','ttensor'}))&&isempty(Ycap)
                    detailsPlot.signal=(Y.U{currentMode})';
                    detailsPlot.num=size((Y.U{currentMode}),2);
                    detailsPlot.stitle=['Source signals - mode ' num2str(currentMode)];
                    detailsPlot.sposition=plotPos1;

                    allPlot.signal=Y.U{currentMode}';
                    allPlot.title=['Source signals - mode ' num2str(currentMode)];
                    allPlot.opos=plotPos2;
                elseif any(strcmpi(class(Y),{'ktensor','ttensor'}))&&(~isempty(Ycap))
                    ROWS=2;
                    detailsPlot(1).signal=(Y.U{currentMode})';
                    detailsPlot(1).num=size((Y.U{currentMode}),2);
                    detailsPlot(1).stitle=['Source signals - mode ' num2str(currentMode)];
                    detailsPlot(1).sposition=plotPos1;

                    allPlot(1).signal=Y.U{currentMode}';
                    allPlot(1).title=['Source signals - mode ' num2str(currentMode)];
                    allPlot(1).opos=plotPos2;

                    detailsPlot(2)=detailsPlot(1);
                    detailsPlot(2).signal=(Ycap.U{currentMode})';
                    detailsPlot(2).num=size((Ycap.U{currentMode}),2);
                    detailsPlot(2).stitle=['Estimated signals - mode ' num2str(currentMode)];
                    detailsPlot(2).sposition=plotPos3;
                    detailsPlot(2).tag='tdvfigsdp2';

                    allPlot(2)=allPlot(1);
                    allPlot(2).signal=Ycap.U{currentMode}';
                    allPlot(2).title=['Estimated signals - mode ' num2str(currentMode)];
                    allPlot(2).opos=plotPos4;
                    allPlot(2).tag='tdvfigsap2';
%                     allPlot(2).ylabel=horzcat('{\bf{a}}^{(',num2str(currentMode),')}_j');
                else
                    errordlg('ttensors are required.','Input error','modal');
                    return;
                end
                for r=1:ROWS
                    detailsPlot(r).sunits=plotUnits;
                    icalab_srplot(detailsPlot(r));
                    h=findobj('tag',horzcat('tdvfigsap',num2str(r)));
                    delete(h);
                    allPlot(r).units=plotUnits;
                    mplotmatrix([],allPlot(r));
                end
            case 'func'
                for r=1:ROWS
                    h=findobj('tag',horzcat('tdvfigsap',num2str(r)));
                    delete(h);
                    allPlot(r).plotfunc=currentPlotFunc;
                    allPlot(r).units=plotUnits;
                    mplotmatrix([],allPlot(r));
                end
        end
    end

    function show3Dtensor
        if NumOfMode<3
            return;
        end
        ROWS=1;      

        if ~isempty(Ycap)&&(~strcmpi(class(Y),'ttensor'))
            core3Dplot.signal=Ycap.core;
            core3Dplot.opos=plotPos1;
            core3Dplot.tag='tdvfigsGcap3D';
            core3Dplot.title='Estimated core tensor';

            if drawSlice
                G3=double(permute(Ycap.core,[Xdims Ydims Zdims]));
                Gsz=size(Ycap.core);
                G3=reshape(G3,prod(Gsz(Xdims)),prod(Gsz(Ydims)),prod(Gsz(Zdims)));
                coreslice.title='Estimated core tensor';
                coreslice.signal=G3;
                coreslice.opos=plotPos2;
                coreslice.tag='tdvfigsGcapslice';
            end
        elseif strcmpi(class(Y),'ttensor')&&isempty(Ycap)
            core3Dplot.signal=Y.core;
            core3Dplot.title='Core tensor';
            core3Dplot.opos=plotPos1;
            core3Dplot.tag='tdvfigsG3D';

            if drawSlice
                G3=double(permute(Y.core,[Xdims Ydims Zdims]));
                Gsz=size(Y.core);
                G3=reshape(G3,prod(Gsz(Xdims)),prod(Gsz(Ydims)),prod(Gsz(Zdims)));
                coreslice.title='Core tensor';
                coreslice.signal=G3;
                coreslice.opos=plotPos2;
                coreslice.tag='tdvfigsGslice';
            end
        elseif strcmpi(class(Y),'ttensor')&&strcmpi(class(Ycap),'ttensor')
            ROWS=2;
            core3Dplot(1).signal=Y.core;
            core3Dplot(1).title='Core tensor';
            core3Dplot(1).opos=plotPos1;
            core3Dplot(1).tag='tdvfigsG3D';

            core3Dplot(2)=core3Dplot(1);
            core3Dplot(2).signal=Ycap.core;
            core3Dplot(2).title='Estimated core tensor';
            core3Dplot(2).opos=plotPos3;
            core3Dplot(2).tag='tdvfigsGcap3D';

            if drawSlice
                G3=double(permute(Y.core,[Xdims Ydims Zdims]));
                Gsz=size(Y.core);
                G3=reshape(G3,prod(Gsz(Xdims)),prod(Gsz(Ydims)),prod(Gsz(Zdims)));
                coreslice(1).title='Core tensor';
                coreslice(1).signal=G3;
                coreslice(1).opos=plotPos2;
                coreslice(1).tag='tdvfigsGslice';

                coreslice(2)=coreslice(1);
                G3=double(permute(Ycap.core,[Xdims Ydims Zdims]));
                Gsz=size(Ycap.core);
                G3=reshape(G3,prod(Gsz(Xdims)),prod(Gsz(Ydims)),prod(Gsz(Zdims)));
                coreslice(2).title='Estimated core tensor';
                coreslice(2).signal=G3;
                coreslice(2).opos=plotPos4;
                coreslice(2).tag='tdvfigsGcapslice';
            end
        else
            errordlg('ttensors are required.','Input error','modal');
            return;
        end
        for r=1:ROWS
            if orderChanged
                core3Dplot(r).units=plotUnits;
                h=findobj('tag','tdvfigsGcap3D');
                if ~isempty(h) 
                    delete(h); 
                end
                h=visual3D({Xdims,Ydims,Zdims},core3Dplot(r));
            end
        end
        orderChanged=false;
        if drawSlice
            h=findobj('tag','tdvfigsGcapslice');
                if ~isempty(h) 
                    delete(h); 
                end
            for r=1:ROWS
                coreslice(r).imagefunc=currentImgFunc;
                coreslice(r).units=plotUnits;
                imagematrix(coreslice(r));
            end
        else
            h=findobj('tag','tdvfigsGslice','-or','tag','tdvfigsGcapslice');
            delete(h);
        end
    end
    function showLambda
%         ROWS=1;
        hbar=CPrank1show(Ycap);
        figure(hfigVMenu);

    end

end %% end of CPTuckerVisualization