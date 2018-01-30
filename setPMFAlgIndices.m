function PMFalgIDs=setPMFAlgIndices
global NumOfMode PMFActParas
global XSpread YSpread backGroundColor lFontSize defaultFontName sFontSize;

if NumOfMode<1
    errordlg('Load a tensor data first.','Input error');
    return;
end
    
PMFActParas=[];
PMFalgIDs=ones(1,NumOfMode);

bkcolor=[0.8314    0.8157    0.7843];
nBorderWidth=.01;
dnBorderWidth=nBorderWidth*2;

PMFActParas=[];

    hfigAlgID=figure('position',[0 0 max(240*XSpread,NumOfMode*150*XSpread) 250*YSpread],'Units','Characters','menu','none','toolbar','none',...
        'Name','Select BSS algorithm','numbertitle','on','color',backGroundColor,'windowstyle','modal','visible','off','CloseRequestFcn',@exitBSSAlgID);
    uicontrol('parent',hfigAlgID,'units','normalized','position',[nBorderWidth 0.8 1-dnBorderWidth 0.15],'style','text',...
        'string','Set BSS algorithms for each mode:','backgroundcolor',bkcolor,'fontsize',lFontSize,'fontname',defaultFontName);
    w=(1-dnBorderWidth)/NumOfMode;wp=w*0.8;
    for n=1:NumOfMode
        uicontrol('parent',hfigAlgID,'units','normalized','position',[nBorderWidth+(n-1)*w 0.65 wp 0.15],'style','text',...
            'string',horzcat('Mode ',num2str(n)),'backgroundcolor',bkcolor,'fontsize',lFontSize,'fontname',defaultFontName,'horizontalalignment','left');
        uicontrol('parent',hfigAlgID,'units','normalized','position',[nBorderWidth+(n-1)*w 0.55 wp 0.15],'style','popupmenu',...
            'string',defstr('penalized matrix factorization'),'fontsize',lFontSize,'fontname',defaultFontName,'horizontalalignment','left','value',1,...
            'userdata',n,'callback',@selPMFalg,'createFcn',@selPMFalg);
    end
    uicontrol('parent',hfigAlgID,'units','normalized','position',[nBorderWidth 0.45 1-dnBorderWidth 0.1],'style','text',...
            'string','*[Free]: do not run BSS on this mode.','backgroundcolor',bkcolor,'fontsize',sFontSize+1,'fontname',defaultFontName,'horizontalalignment','left');
    uicontrol('parent',hfigAlgID,'units','normalized','position',[0.3 0.3 0.4 0.13],'style','pushbutton',...
            'string','Algorithm options','fontsize',lFontSize,'fontname',defaultFontName,'callback',@setPMFalgOpts);    
    uicontrol('parent',hfigAlgID,'units','normalized','position',[0.4 0.05 0.2 0.15],'style','pushbutton',...
            'string','Confirm','fontsize',lFontSize,'fontname',defaultFontName,'horizontalalignment','left','value',1,...
            'callback',@exitBSSAlgID);
    movegui(hfigAlgID,'center');
    set(hfigAlgID,'visible','on');
    waitfor(hfigAlgID);
            
                
        function selPMFalg(hObject,eventdata,source)
            PMFalgIDs(get(hObject,'userdata'))=get(hObject,'value');
        end
        function setPMFalgOpts(hObject,eventdata,source)
            flag=PMFalgIDs>1;
            if sum(flag)>0
                %% call guiSetOpts                  
                  [PMFparaList PMFalgs PMFparaTList]=PMFalgInit();
                  algPMFalgIDs=PMFalgIDs(flag)-1;
                  opts=guiSetOpts(PMFalgs(algPMFalgIDs),PMFparaTList(algPMFalgIDs),PMFparaList(algPMFalgIDs));
                  pos=1;
                  PMFActParas=cell(1,numel(algPMFalgIDs));
                  for cidx=1:numel(algPMFalgIDs)
                      PMFActParas{cidx}=opts{pos};
                      pos=pos+1;
                  end
            end
        end
    function exitBSSAlgID(hObject,eventdata,source)
        PMFalgIDs=reshape(PMFalgIDs,1,[]);

        %% set default options for PMF here if it is empty
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

        delete(hfigAlgID);
    end



end