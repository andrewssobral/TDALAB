function genCoreTensor(hObject,eventdata,handles)
global Y NumOfComp NumOfMode tdalabStatus;
global ScreenWidth ScreenHeight XSpread YSpread lFontSize defaultFontName backGroundColor;
% hObject=findobj('tag','pmG');
Gmode=get(hObject,'value');
if isempty(tdalabStatus.inputType)
    errordlg('No tensor/ktensor/ttensor found.','No input','modal');
    return;
elseif strcmpi(tdalabStatus.inputType,'tensor')
    errordlg('Only valid for ttensor/ktensor.','Error','modal');
    return;
end
GGenMode=defstr('gencore');

hTType=findobj('tag','pmInputFormat');
TType=get(hTType,'string');
if strcmpi(TType{get(hTType,'value')},'ttensor')
% if strcmp(tdalabStatus.model,'Tucker')
    if Gmode==1
        return;
    end
    sp=0;
    H=13*YSpread;W=70*XSpread;
    hslider=[];
    
    if numel(NumOfComp)==1
        NumOfComp=NumOfComp(ones(1,NumOfMode));
    end

    hmainSPL=figure('Units','Characters','toolbar','none','menu','none','name',horzcat('Generate core tensor [',GGenMode{Gmode},' ]'),...
        'position',[ScreenWidth/2-W/2 ScreenHeight/2 W H],'numbertitle','off','resize','on',...
        'WindowStyle','modal','color',backGroundColor);
    bkcolor=get(hmainSPL,'color');
    nBorderWidth=0.05;

    htxtTitle=uicontrol('parent',hmainSPL,'Units','normalized','position',[nBorderWidth 0.75 1-2*nBorderWidth 0.15],...
        'style','text','string','Sparsity of core tensor [G]','backgroundcolor',bkcolor,...
        'fontsize',lFontSize);

    hslider=uicontrol('parent',hmainSPL,'Units','normalized','position',[nBorderWidth 0.6 1-2*nBorderWidth 0.13],...
        'style','slider','Max',1,'Min',0,'foregroundcolor','blue','sliderstep',[0.01 0.1],...
        'fontsize',lFontSize,'tag','sliderspa','callback',@UpdateSlider,...
        'backgroundcolor',[0.3 0.7 0.3],'value',sp);

    htxtValue=uicontrol('parent',hmainSPL,'Units','normalized','position',[0.75 0.45 0.25-nBorderWidth*2 0.15],...
        'style','text','string',strcat(num2str(round(100*sp)),'%'),'foregroundcolor','blue','backgroundcolor',bkcolor,...
        'fontsize',lFontSize,'HorizontalAlignment','right');
    
    hckParaTuck=uicontrol('parent',hmainSPL,'units','normalized','position',[0.6 0.3 0.39 0.15],...
        'style','checkbox','string','ParaTuck2','backgroundcolor',bkcolor,'fontsize',lFontSize,'horizontalalignment','right',...
        'fontname',defaultFontName,'enable','off','Max',true,'Min',false);

    hpbOK=uicontrol('parent',hmainSPL,'Units','normalized','position',[0.15 0.05 0.25 0.17],...
        'style','pushbutton','string','OK','callback',@pbOK,...
        'fontsize',lFontSize);
    hpbCancel=uicontrol('parent',hmainSPL,'Units','normalized','position',[0.6 0.05 0.25 0.17],...
        'style','pushbutton','string','Cancel','callback',@pbCancel,...
        'fontsize',lFontSize);
    
    if NumOfComp(1)==NumOfComp(2)&&NumOfMode==3
        set(hckParaTuck,'enable','on','value',false);
    end
else
    switch Gmode
        case {1,2}
            lambda=rand(NumOfComp(1),1);
        case 3
            lambda=randn(NumOfComp(1),1);
        case 4
            lambda=abs(randn(NumOfComp(1),1));
    end
    Y=ktensor(lambda,Y.U);
    cleartemp;
end

        function UpdateSlider(hObject, eventdata, handles)
            sp=get(hObject,'value');
            str=[num2str(round(sp*100)),'%'];
            set(htxtValue,'string',str);
        end
        function pbOK(hObject, eventdata, handles)
            sp=get(hslider,'value'); 
            switch Gmode
                case 2
                    G=rand(NumOfComp);
                    G=G.*(G>sp);
                case 3
                    G=randn(NumOfComp);
                    G=G.*(rand(NumOfComp)>sp);
                case 4
                    G=abs(randn(NumOfComp));
                    G=G.*(rand(NumOfComp)>sp);
            end
            
            if get(hckParaTuck,'value')
                ind=zeros(NumOfComp(3)*NumOfComp(1),3);
                z=(1:NumOfComp(1))';
                ind(:,1:2)=repmat(reshape(z(:,ones(1,NumOfComp(3)))',size(ind,1),1),1,2);
                ind(:,3)=repmat((1:NumOfComp(3))',NumOfComp(1),1);
                G=double(G).*double(sptensor(ind,1));
            end
            
            Y=ttensor(tensor(G),Y.U);
            cleartemp;
            delete(hmainSPL);
        end
        function pbCancel(hObject, eventdata, handles)
            delete(hmainSPL);
        end



end

