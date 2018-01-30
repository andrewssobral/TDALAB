function visualize
global tdalabStatus XSpread YSpread algs;
global Y Ycap MCSIRs NumOfMode backGroundColor lFontSize defaultFontName;
global MCfit MCelapsedTime;

if ~tdalabStatus.decomposed
    errordlg('Please run tensor decomposition algorithm first.','No input','modal');
    return;
end

plotFuncstr={'plot3','plot','waterfall','mplot','stem3','ribbon'};


ctrlHeight=2.2*YSpread;
txtHeight=1.6*YSpread;
ctrlWidth=20*XSpread;

if tdalabStatus.advEvaluation
    MCVisualize;
    return;
end

tdalab('hide');

switch tdalabStatus.model
    case {'CP','Tucker','PMF'}
        CPTuckerVisualize;
    case 'BCD'
        BCDVisualize;
    otherwise
        errordlg('Unsuported tensor decoposition model.','Model error','modal');
        return;
end

%     function parafac2Visualize
%         hpmFunc=[];currentPlotFunc='plot3';
%         h=findobj('tag','mtdvfigMenuPara2');
%         if ~isempty(h)
%             figure(h);
%             h=findobj('-regexp','tag','\<tdvfig.*');
%             set(h,'visible','on');
%             return;
%         end
%         delete(findobj('-regexp','tag','\<tdvfig.*'));
%         switch NumOfMode
%             case 3
%                 hfmparafac2=figure('units','normalized','outerposition',[0. 0. 0.1 0.15],'name','Parafac2','numbertitle','off',...
%                     'toolbar','none','menu','none','color',backGroundColor,'tag','mtdvfigMenuPara2');
%                 movegui(hfmparafac2,'northwest');
%                 uicontrol('parent',hfmparafac2,'units','normalized','position',[0.01 0.7 0.98 0.2],'style','text',...
%                     'backgroundcolor',backGroundColor,'Fontsize',lFontSize,'fontname',defaultFontName,'string','Draw function');
%                 hpmFunc=uicontrol('parent',hfmparafac2,'units','normalized','position',[0.01 0.5 0.98 0.2],'style','popupmenu',...
%                     'Fontsize',lFontSize,'fontname',defaultFontName,'string',plotFuncstr,'callback',@parafac2Func);
%                 uicontrol('parent',hfmparafac2,'units','normalized','position',[0.01 0.01 0.98 0.3],'style','pushbutton',...
%                     'Fontsize',lFontSize,'fontname',defaultFontName,'string','TDALAB','callback','tdalab;');
%                 drawParafac2;
%             otherwise
%         end
%         function parafac2Func(hObject,eventdata)
%             currentPlotFunc=plotFuncstr{get(hObject,'value')};
%             drawParafac2;
%         end
% 
%         function drawParafac2
%             delete(findobj('tag','tdvfigs_sub'));
%             set(hfmparafac2,'units','normalized');
%             pos=get(hfmparafac2,'outerposition');
%             left=pos(1)+pos(3)+0.01;
%             top=pos(2)+pos(4);
%             if strcmpi(tdalabStatus.inputType,'tensor')
%                 height=0.6;width=(1-left);
%                 hm=figure('units','normalized','outerposition',[left top-height width height],'name','Parafac2','numbertitle','off',...
%                     'color',backGroundColor,'tag','tdvfigs_sub');
%                 movegui(hm,'center');
%                 info(1)=struct('signal',Ycap.A','title','First mode','plotfunc',currentPlotFunc);
%                 info(2)=struct('signal',Ycap.C','title','Third mode','plotfunc',currentPlotFunc);
%                 info(3)=struct('signal',(Ycap.P{1}*Ycap.H)','title','Second mode (first k-slab)','plotfunc',currentPlotFunc);
%                 info(4)=struct('signal',(Ycap.P{2}*Ycap.H)','title','Second mode (second k-slab)','plotfunc',currentPlotFunc);
%                 for n=1:4
%                     hax=subplot(2,2,n);
%                     mplotmatrix(hax,info(n));
%                 end
%             else
%                 width=1-left;height=0.48;
%                 figure('units','normalized','outerposition',[left top-height width height],'name','Parafac2 - estimates','numbertitle','off',...
%                     'color',backGroundColor,'tag','tdvfigs_sub');
%                 info(1)=struct('signal',Ycap.A','title','First mode','plotfunc',currentPlotFunc);
%                 info(2)=struct('signal',Ycap.C','title','Third mode','plotfunc',currentPlotFunc);
%                 info(3)=struct('signal',(Ycap.P{1}*Ycap.H)','title','Second mode (first k-slab)','plotfunc',currentPlotFunc);
%                 info(4)=struct('signal',(Ycap.P{2}*Ycap.H)','title','Second mode (second k-slab)','plotfunc',currentPlotFunc);
%                 for n=1:4
%                     hax=subplot(2,2,n);
%                     mplotmatrix(hax,info(n));
%                 end
%                 figure('units','normalized','outerposition',[left top-height*2 width height],'name','Parafac2 - sources','numbertitle','off',...
%                     'color',backGroundColor,'tag','tdvfigs_sub');
%                 info(1)=struct('signal',Y.U{1}','title','First mode','plotfunc',currentPlotFunc);
%                 info(2)=struct('signal',Y.U{3}','title','Third mode','plotfunc',currentPlotFunc);
%                 info(3)=struct('signal',(Y.U{2})','title','Second mode','plotfunc',currentPlotFunc);
%                 for n=1:3
%                     hax=subplot(2,2,n);
%                     mplotmatrix(hax,info(n));
%                 end
%                 if ~isempty(MCSIRs)
%                     showPIs;
%                 end
%             end
%         end
%     end

    function BCDVisualize        
        hfigsBCD=[];currentPlotFunc='plot3';
        currentBlock=1;
        blockChanged=true;
        h=findobj('tag','mtdvfigMenuBCD');
        if ~isempty(h)
            figure(h);
            h=findobj('-regexp','tag','\<tdvfig.*');
            set(h,'visible','on');
            return;
        end
        delete(findobj('-regexp','tag','\<tdvfig.*'));
        
        if iscell(Ycap.A)
            NBlocks=size(Ycap.A,2);
        else
            NBlocks=1;
        end
        hfmBCD=figure('units','characters','position',[0 0 ctrlWidth+2*XSpread 12.5*YSpread],'name','BCD','numbertitle','off',...
                    'toolbar','none','menu','none','color',backGroundColor,'tag','mtdvfigMenuPara2');
        movegui(hfmBCD,'northwest');
        uicontrol('parent',hfmBCD,'units','characters','position',[XSpread 10.2*YSpread  ctrlWidth txtHeight],'style','text',...
            'backgroundcolor',backGroundColor,'fontsize',lFontSize,'fontname',defaultFontName,'string','Block');
        uicontrol('parent',hfmBCD,'units','characters','position',[XSpread 8*YSpread  ctrlWidth ctrlHeight],'style','popupmenu',...
            'fontsize',lFontSize,'fontname',defaultFontName,'string',1:NBlocks,'callback',@changeBlock);
        uicontrol('parent',hfmBCD,'units','characters','position',[XSpread 5.7*YSpread  ctrlWidth txtHeight],'style','text',...
            'backgroundcolor',backGroundColor,'fontsize',lFontSize,'fontname',defaultFontName,'string','Plot function');
        uicontrol('parent',hfmBCD,'units','characters','position',[XSpread 3*YSpread  ctrlWidth ctrlHeight],'style','popupmenu',...
            'fontsize',lFontSize,'fontname',defaultFontName,'string',plotFuncstr,'callback',@BCDchangePlot);
        uicontrol('parent',hfmBCD,'units','characters','position',[XSpread 0.1*YSpread  ctrlWidth ctrlHeight],'style','pushbutton',...
            'fontsize',lFontSize,'fontname',defaultFontName,'string','TDALAB','callback','tdalab;');
        BCDVupdate;
        
        function changeBlock(hObject,eventdata)
            currentBlock=get(hObject,'value');
            blockChanged=true;
            BCDVupdate;
        end
        function BCDchangePlot(hObject,eventdata)
            currentPlotFunc=plotFuncstr{get(hObject,'value')};
            blockChanged=false;
            BCDVupdate;
        end
        function BCDVupdate
            oldUnits=get(0,'units');
            set(0,'units','inch');
            WinSize=get(0,'monitorposition');
            WinSize=WinSize(1,:);
            set(0,'units',oldUnits);
            set(hfmBCD,'units','inch');            
            pos=get(hfmBCD,'outerposition');
            width=min(WinSize(3)-sum(pos([1 3]))-0.01,8);
            height=min(WinSize(4),width);
            plotPos=[sum(pos([1 3]))+0.01 sum(pos([2,4]))-height width height];
%                     plotPos=[sum(pos([1 3]))+0.01 sum(pos([2,4]))-0.8 0.98-sum(pos([1 3])) 0.8];
                    h=findobj('tag','tdvfigsBCD');
                    if isempty(h)
                        figure('units','inch','outerposition',plotPos,'name',horzcat('BCD-',Ycap.type),'numbertitle','off',...
                            'color',backGroundColor,'tag','tdvfigsBCD');
                    else
                        figure(h);
                    end
            switch Ycap.type
                case {'bcdLMN','bcdLrMrNr'}
                    info(1)=struct('signal',Ycap.A{currentBlock}','title',horzcat('{\bfA}^{(',num2str(currentBlock),')}'),'tag','tdvfigsBDCsub','plotfunc',currentPlotFunc);
                    info(2)=struct('signal',Ycap.B{currentBlock}','title',horzcat('{\bfB}^{(',num2str(currentBlock),')}'),'tag','tdvfigsBDCsub','plotfunc',currentPlotFunc);
                    info(3)=struct('signal',Ycap.C{currentBlock}','title',horzcat('{\bfC}^{(',num2str(currentBlock),')}'),'tag','tdvfigsBDCCore','plotfunc',currentPlotFunc);    
                    for i=1:3haxes=subplot(2,2,i);
                        mplotmatrix(haxes,info(i));
                    end
                    if blockChanged
                        coreinfo=struct('signal',Ycap.D{currentBlock},'title',horzcat('Core tensor: {\itD}_',num2str(currentBlock)));
                        haxes=subplot(2,2,4);
                        coreinfo.opos=get(haxes,'outerposition');
                        coreinfo.tag='tdvfigsBDCsub';
                        visual3D({[1],[2],[3]},coreinfo,haxes);
                    end
                case {'bcdLL1','LrLr1'}
                        info(1)=struct('signal',Ycap.A','title','\bfA','tag','tdvfigsBDCsub','plotfunc',currentPlotFunc);
                        info(2)=struct('signal',Ycap.B','title','\bfB','tag','tdvfigsBDCsub','plotfunc',currentPlotFunc);
                        info(3)=struct('signal',Ycap.C','title','\bfC','tag','tdvfigsBDCCore','plotfunc',currentPlotFunc);
                    for i=1:3
                        haxes=subplot(2,2,i);
                        mplotmatrix(haxes,info(i));
                    end
                case {'bcdLM','bcdLrMr'}
                        info(1)=struct('signal',Ycap.A{currentBlock}','title',horzcat('{\bfA}^{(',num2str(currentBlock),')}'),'tag','tdvfigsBDCsub','plotfunc',currentPlotFunc);
                        info(2)=struct('signal',Ycap.B{currentBlock}','title',horzcat('{\bfB}^{(',num2str(currentBlock),')}'),'tag','tdvfigsBDCsub','plotfunc',currentPlotFunc);
                    for i=1:2
                        haxes=subplot(2,2,i);
                        mplotmatrix(haxes,info(i));
                    end
                    if blockChanged
                        coreinfo=struct('signal',Ycap.C{currentBlock},'title',horzcat('Core tensor: {\itC}_',num2str(currentBlock)));
                        haxes=subplot(2,2,3);
                        visual3D({[1],[2],[3]},coreinfo,haxes);
                    end
                otherwise
                    errordlg('Unexpected return type.','Unknown error','modal');
            end
        end
    end

end