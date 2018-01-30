function updateUI
global tdalabStatus Y NumOfComp NumOfMode haxesTensorInfor defaultFontName;
inputFormatstr=defstr('input');
TDModelstr=defstr('TDModel');
showTensorInfor;
if tdalabStatus.fullTensor
    h=findobj('tag','txtMode','-or','tag','txtCoreTensor',...
        '-or','tag','ppMode','-or','tag','pmG','-or','tag','pbModi',...
        '-or','tag','txtInputFormat','-or','tag','pmInputFormat');
    set(h,'enable','off');
    
    h=findobj('tag','txtNoise','-or','tag','pmNoise');
    set(h,'enable','on');
    
    htxt=findobj(h,'tag','txtNoise');
    if isinf(tdalabStatus.noiseSNR)
        noisestr=horzcat('Add noise [Current: noise free]');
    else
        noisestr=horzcat('Add noise [',tdalabStatus.noiseType, ', SNR=',num2str(tdalabStatus.noiseSNR),' dB]');
    end
    set(htxt,'string',noisestr);
%     htype=findobj(h,'tag','pmNoise');
%     idx=find(strcmpi(get(htype,'string'),tdalabStatus.noiseType));
%     if isempty(idx)
%         set(htype,'value',1);
%     else
%         set(htype,'value',idx);
%     end
    
    
    if strcmp(tdalabStatus.inputType,'tensor')
        h=findobj('tag','ckToTensor');
        set(h,'value',true,'enable','off');
    else
        set(h,'value',true,'enable','on');
    end
    hmodel=findobj('tag','pmTDModel');
    if isempty(tdalabStatus.model)
        tdalabStatus.model=TDModelstr{get(hmodel,'value')};
    else        
        n=find(strcmpi(TDModelstr,tdalabStatus.model),1);
        set(hmodel,'value',n);
    end
else
    h=findobj('tag','txtMode','-or','tag','txtCoreTensor',...
        '-or','tag','ppMode','-or','tag','pmG','-or','tag','pbModi',...
         '-or','tag','txtInputFormat','-or','tag','pmInputFormat');
    set(h,'enable','on');
    hInputFormat=findobj(h,'tag','pmInputFormat');
    if ~isempty(tdalabStatus.inputType)
        v=find(strcmpi(inputFormatstr,tdalabStatus.inputType),1);
        if isempty(v)
            errordlg('Unexpected input type.','Input error','modal');
            return;
        else
            set(hInputFormat,'value',v);
        end
    end
    hG=findobj(h,'tag','txtCoreTensor');
    if strcmpi(tdalabStatus.inputType,'ktensor')
        set(hG,'string','lambda:');
        h=findobj('tag','gbExtraConstraints');
        set(allchild(h),'enable','off');
    elseif strcmpi(tdalabStatus.inputType,'ttensor')
        set(hG,'string','Core tensor:');
        hmodel=findobj('tag','pmTDModel');
        if strcmpi(TDModelstr{get(hmodel,'value')},'Tucker')
            h=findobj('tag','gbExtraConstraints');
            set(allchild(h),'enable','on');
        end
    end
    h=findobj('tag','ppMode');set(h,'string',{1:NumOfMode});
    set(findobj('tag','txtNoise','-or','tag','pmNoise'),'enable','off');
    h=findobj('tag','ckToTensor');
    set(h,'value',false,'enable','on');
end

if isempty(tdalabStatus.inputType)
    h=findobj('tag','pbRunAlg');
    set(h,'enable','off');
else
    h=findobj('tag','pbRunAlg');
    set(h,'enable','on');
end

h=findobj('tag','gbExtraConstraints');
if strcmpi(tdalabStatus.model,'Tucker')
    set(allchild(h),'enable','on');
    set(h,'SelectedObject',[]);
else
    set(allchild(h),'enable','off');
end

if ~isempty(tdalabStatus.inputType)&&any(strcmpi({'Tucker','CP'},tdalabStatus.model))
    set(findobj('tag','pbMCRunOptions','-or','tag','txtAdvEvaluation'),'enable','on');
else
    set(findobj('tag','pbMCRunOptions','-or','tag','pbMCRun','-or','tag','txtAdvEvaluation'),'enable','off');
end

h=findobj('tag','cbNN');
set(h,'value',tdalabStatus.nonnegativity);

if tdalabStatus.decomposed
    set(allchild(findobj('tag','pnlOutputAnalysis')),'enable','on');
else
    set(allchild(findobj('tag','pnlOutputAnalysis')),'enable','off');
end
set(findobj('tag','pbgotoCCwin'),'enable','on');

updateAlgList;


    function showTensorInfor
        h=findobj('tag','txtNoTensorInfor');
        if isempty(tdalabStatus.inputType)
            set(h,'visible','on');
            return;
        end
        set(h,'visible','off');
        
        axes(haxesTensorInfor);
        cla(haxesTensorInfor);
        axis('off');
        xlim([0 1]);ylim([0 1]);
        szY=size(Y);
        if ~isempty(szY)
            szYstr=num2str(szY(1));
            for n=2:NumOfMode
                szYstr=horzcat(szYstr,'$\times$',num2str(szY(n)));
            end
            
            %\underline{\bf{Y}}','interpreter','latex'
            szYstr=horzcat('Dim. of {\underline{\bf{Y}}}: ',szYstr);
            if strcmp(tdalabStatus.model,'Tucker')
                if (NumOfComp(1)>=1)&&all(NumOfComp>0)
                    szGstr=num2str(NumOfComp(1));
                    for k=2:length(NumOfComp)
                        szGstr=horzcat(szGstr,'$\times$',num2str(NumOfComp(k)));
                    end
                else
                    szGstr='Unknown.';
                end
                szGstr=horzcat('Dim. of {\underline{\bf{G}}}: ',szGstr);
            elseif strcmp(tdalabStatus.model,'CP')
                if all(NumOfComp==0)
                    szGstr=horzcat('Rank [NumOfComp]: Unknown.');
                else
                    szGstr=horzcat('Rank [NumOfComp]: ',num2str(NumOfComp(1)));
                end
            else
                szGstr='Core tensor: Unknown.';
            end
            text(0.005,0.85,horzcat('Type of {\underline{\bf{Y}}}: ',num2str(NumOfMode),'-way ',class(Y)),'fontname',defaultFontName,'interpreter','latex');
            text(0.005,0.5,szYstr,'interpreter','latex');
            text(0.005,0.15,szGstr,'interpreter','latex');
        else
            text(0.005,0.5,'No Input.');
        end

    end
end