function multiRunOptions
global algs tdalabStatus MCTimes;
global XSpread YSpread backGroundColor lFontSize defaultFontName ;

Width=120*XSpread;
Height=30*YSpread;
fmainMultiRun=figure('Units','Characters','Resize','on','toolbar','none','menu','none',...
    'Position',[0,0,Width,Height],'Windowstyle','modal','numbertitle','off',...
    'Name','Settings ...','color',backGroundColor);
movegui(fmainMultiRun,'center');

validNames={algs(tdalabStatus.validAlgs).details};
Nvalg=numel(tdalabStatus.validAlgs);
leftList=(1:Nvalg)';
rightList=[];

hpnlAlgs=uipanel('parent',fmainMultiRun,'units','normalized','position',[0.01 0.25 0.98 0.72],'title','','backgroundcolor',backGroundColor,...
    'borderwidth',1);
uicontrol('parent',hpnlAlgs,'units','normalized','position',[0.01 0.85 0.449 0.12],'style','text','string','Valid algorithms',...
    'backgroundcolor',backGroundColor,'fontsize',lFontSize,'fontname',defaultFontName);
uicontrol('parent',hpnlAlgs,'units','normalized','position',[0.551 0.85 0.449 0.12],'style','text','string','Selected algorithms',...
    'backgroundcolor',backGroundColor,'fontsize',lFontSize,'fontname',defaultFontName);
hleftList=uicontrol('parent',hpnlAlgs,'units','normalized','position',[0.01 0.01 0.448 0.84],'style','list','string',validNames(leftList),...
    'backgroundcolor','white','fontsize',lFontSize,'fontname',defaultFontName,'max',1,'min',1,'value',1);
hrightList=uicontrol('parent',hpnlAlgs,'units','normalized','position',[0.55 0.01 0.44 0.84],'style','list','string','',...
    'backgroundcolor','white','fontsize',lFontSize,'fontname',defaultFontName,'max',1,'min',1,'value',0);
haddAlg=uicontrol('parent',hpnlAlgs,'units','normalized','position',[0.46 0.55 0.09 0.1],'style','pushbutton','string','>>',...
    'fontsize',lFontSize,'fontname',defaultFontName,'fontweight','bold','callback',@addOneAlg);
hrmAlg=uicontrol('parent',hpnlAlgs,'units','normalized','position',[0.46 0.35 0.09 0.1],'style','pushbutton','string','<<',...
    'fontsize',lFontSize,'fontname',defaultFontName,'fontweight','normal','callback',@removeOneAlg);

uicontrol('parent',fmainMultiRun,'units','normalized','position',[0.01,0.13 0.2 0.05],'style','text','string','Repeated times:',...
    'backgroundcolor',backGroundColor,'fontsize',lFontSize,'fontname',defaultFontName,'horizontalalign','left');
hMCTimes=uicontrol('parent',fmainMultiRun,'units','normalized','position',[0.21,0.13 0.3 0.07],'style','edit','string','10',...
    'backgroundcolor','white','fontsize',lFontSize,'fontname',defaultFontName,'horizontalalign','left','callback',@checkMCTimes);
hOK=uicontrol('parent',fmainMultiRun,'units','normalized','position',[0.3 .01 0.15 0.07],'style','pushbutton','string','OK',...
    'fontsize',lFontSize,'fontname',defaultFontName,'horizontalalign','left','callback',@OKMultiRun,'enable','off');
hCancel=uicontrol('parent',fmainMultiRun,'units','normalized','position',[0.55 .01 0.15 0.07],'style','pushbutton','string','Cancel',...
    'fontsize',lFontSize,'fontname',defaultFontName,'horizontalalign','left','callback',@cancelMultiRun);

    function addOneAlg(hObject,eventdata)
        sel=get(hleftList,'value');
        leftList(sel)=[];
        set(hleftList,'string',validNames(leftList));
        if isempty(leftList)
            set(haddAlg,'enable','off');
        else
            set(hleftList,'value',1);
        end
        rightList=setdiff(1:Nvalg,leftList)';
        set(hrightList,'string',validNames(rightList));
        set(hrightList,'value',1);
        set([hrmAlg hOK],'enable','on');
    end
    function removeOneAlg(hObject,eventdata)
        sel=get(hrightList,'value');
        rightList(sel)=[];
        set(hrightList,'string',validNames(rightList));
        if isempty(rightList)
            set([hrmAlg hOK],'enable','off');
        else
            set(hrightList,'value',1);
        end
        leftList=setdiff(1:Nvalg,rightList);
        set(hleftList,'string',validNames(leftList));
        set(hleftList,'value',1);
        set(haddAlg,'enable','on');
    end

    function OKMultiRun(hObject,eventdata)
        MCTimes=str2num(get(hMCTimes,'string'));
        button = questdlg('Press Yes to set parameters for each algorithm or No to use default settings.','Please confirm');
        switch button
            case 'Yes'
                tdalabStatus.algIndex=tdalabStatus.validAlgs(rightList);
                delete(fmainMultiRun);
                
                setOpts(tdalabStatus.algIndex);
%                 for i=tdalabStatus.algIndex
%                     h=setOpts(i);
%                     uiwait(h);
%                 end
                
                set(findobj('tag','pbMCRun'),'enable','on');
                tdalabStatus.advEvaluation=true;
                tdalab;
            case 'No'
                tdalabStatus.algIndex=tdalabStatus.validAlgs(rightList);
                delete(fmainMultiRun);
                
                set(findobj('tag','pbMCRun'),'enable','on');
                tdalabStatus.advEvaluation=true;
                tdalab;
            otherwise
                figure(fmainMultiRun);
        end
    end
    function cancelMultiRun(hObject,eventdata)
        delete(fmainMultiRun);
        tdalab;
    end
    function checkMCTimes(hObject,eventdata)
        n=round(str2double(get(hObject,'string')));
        if n<1
            errordlg('This value must be a positive integer.','Invalid input','modal');
            set(hObject,'string','10');
        else
            set(hObject,'string',num2str(n));
        end
    end
end