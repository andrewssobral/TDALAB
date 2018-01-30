function setNoise( hObject,eventdata,handles)
global ScreenWidth ScreenHeight XSpread YSpread  lFontSize defaultFontName;
global tdalabStatus Y TEMPFILE;
NoiseStr={'no noise';'Uniform';'Gaussian';'Abs. Gaussian'};
noiseType=get(hObject,'string');
noiseType=noiseType{get(hObject,'value')};
defaultNoiseLevel=20;NoiseLevel=20;
if ~isempty(tdalabStatus.noiseSparsity)
    NZRatio=tdalabStatus.noiseSparsity;
else
    NZRatio=0;
end
OldNoiseLevel=tdalabStatus.noiseSNR;
nvalue=get(hObject,'value');
if nvalue==1
    tdalabStatus.noiseType=[];
    tdalabStatus.noiseSNR=inf;
    tdalabStatus.noiseSparsity=0;
    cleartemp;
    updateUI;
    return;
end

nBorderWidth=0.05;
Width=70*XSpread;Height=13*YSpread;
hmainNoise=figure('Units','Characters','position',[ScreenWidth/2-Width/2,ScreenHeight/2-Height/2,Width,Height],...
    'toolbar','none','menu','none','resize','on','Name',horzcat('Add noise: [',NoiseStr{nvalue},']'),'numbertitle','off','WindowStyle','modal');
bkcolor=get(hmainNoise,'color');
htxtNoiseLevel=uicontrol('parent',hmainNoise,'Units','normalized','position',[nBorderWidth 0.76 1-2*nBorderWidth 0.15],...
    'style','text','string',horzcat('Noise level ( ',num2str(min(OldNoiseLevel,defaultNoiseLevel),'%d'),' dB)'),'fontname',defaultFontName,'fontsize',lFontSize,'backgroundcolor',bkcolor);
hsdrNoiseLevel=uicontrol('parent',hmainNoise,'Units','normalized','position',[nBorderWidth 0.65,1-2*nBorderWidth 0.15],...
    'style','slider','Max',60,'Min',0,'sliderstep',[1/60,1/12],'foregroundcolor','blue','fontsize',lFontSize,'tag','sdrNoiseLevel',...
    'callback',@setNoiseLevel,'value',min(OldNoiseLevel,defaultNoiseLevel),'backgroundcolor',[0.3 0.7 0.3]);

htxtSpa=uicontrol('parent',hmainNoise,'Units','normalized','position',[nBorderWidth 0.4 1-2*nBorderWidth 0.15],...
    'style','text','string',horzcat('Number of zeros in noise term ( ',num2str(round(NZRatio*100),'%d%%'),' )'),'fontname',defaultFontName,'fontsize',lFontSize,'backgroundcolor',bkcolor);
hsdrSpRatio=uicontrol('parent',hmainNoise,'Units','normalized','position',[nBorderWidth 0.3 1-2*nBorderWidth 0.15],...
    'style','slider','Max',1,'Min',0,'sliderstep',[0.01,0.05],'foregroundcolor','blue','fontsize',lFontSize,'tag','sdrNoiseLevel',...
    'callback',@setSpRatio,'value',NZRatio,'backgroundcolor',[0.3 0.7 0.3]);

hpbOK=uicontrol('parent',hmainNoise,'Units','normalized','position',[0.2 0.05 0.2 0.15],...
    'style','pushbutton','string','OK','fontname',defaultFontName,'fontsize',lFontSize,'callback',@pbOK);
hpbCancel=uicontrol('parent',hmainNoise,'Units','normalized','position',[0.6 0.05 0.2 0.15],...
    'style','pushbutton','string','Cancel','fontname',defaultFontName,'fontsize',lFontSize,'callback',@pbCancel);

    function setNoiseLevel(hObject,eventdata,handles)
        NoiseLevel=floor(get(hObject,'value'));
        set(htxtNoiseLevel,'string',horzcat('Noise level ( ',num2str(NoiseLevel,'%d'),' dB)'));
    end
    function setSpRatio(hObject,eventdata,handles)
        NZRatio=round(100*get(hObject,'value'))/100;
        set(htxtSpa,'string',horzcat('Number of zeros in noise term ( ',num2str(round(NZRatio*100),'%d%%'),' )'));
    end
% NoiseStr={'no noise';'Uniform';'Gaussian';'Abs. Gaussian'};
    function pbOK(hObject,eventdata,handles)
        delete(hmainNoise);
        commandwindow;
        
        hNoisetxt=findobj('tag','txtNoise');
        tdalab('hide');
        tdalabStatus.noiseType=noiseType;
        tdalabStatus.noiseSNR=NoiseLevel;
        tdalabStatus.noiseSparsity=NZRatio;
        fprintf('[TDALAB] Begin to add noise. This may cost a few minutes depending on the size of input.\n');
        fprintf('[TDALAB] Please wait and do not press any key...\n');
        pause(0.1);
        addNoise;
        fprintf('[TDALAB] Noise is added. SNR: %ddB;\tSparsity: %d%%;\tType: %s.\n',NoiseLevel,tdalabStatus.noiseSparsity*100,noiseType);
        
        
        if isinf(tdalabStatus.noiseSNR)
            noisestr=horzcat('Add noise [Current: noise free]');
        else
            noisestr=horzcat('Add noise [Current: ',tdalabStatus.noiseType,'',num2str(tdalabStatus.noiseSNR),'dB]');
        end
        set(hNoisetxt,'string',noisestr);
        
        tdalab('show');
    end
    function pbCancel(hObject,eventdata,handles)
        delete(hmainNoise);
    end

    function addNoise
        Ydim=size(Y);
        
%         L=prod(Ydim);
%         flag=rand(1,L);
%         flag=flag<NZRatio;
%         switch nvalue
%             case 2 % Uniform
%                 noi=rand(1,L);
%                 noi(flag)=0;
%             case 3 % Gaussian
%                 noi=randn(1,L);
%                 noi(flag)=0;
%             case 4 % abs. Gaussian
%                 noi=abs(randn(1,L));
%                 noi(flag)=0;
%             otherwise
%                 errordlg('Error noise type. Fail to add noise.','Error','modal');
%                 return;
%         end        
        
%         Yfull=double(Y);
%         Yfullvec=reshape(Yfull,1,L);
%         Yfullp=Yfullvec*Yfullvec';
%         noip=noi*noi';        
%         r=sqrt(Yfullp/noip)*10^(-NoiseLevel/20);
%         Yfullvec=Yfullvec+r*noi;
%         clear noi Yfull;
%         % backup
%         Ynoise=tensor(reshape(Yfullvec,Ydim));
        
        
        fprintf('.');pause(0.1);
        
        switch nvalue
            case 2 % Uniform
                noi=rand(Ydim);
            case 3 % Gaussian
                noi=randn(Ydim);
            case 4 % abs. Gaussian
                noi=abs(randn(Ydim));
            otherwise
                errordlg('Error noise type. Fail to add noise.','Error','modal');
                return;
        end

        flag=rand(Ydim);
        flag=flag<NZRatio;
        
        noi(flag)=0;
        clear flag;
        
        fprintf('.');pause(0.1);
        
        noi=tensor(noi);
        noip=norm(noi)^2;        
        Yfull=tensor(Y);
        Yfullp=norm(Yfull)^2;
        
        r=sqrt(Yfullp/noip)*10^(-NoiseLevel/20);
        
        Ynoise=Yfull+r.*noi;
        clear noi Yfull;
        
        fprintf('.\n');pause(0.1);
        
        
        SNR=tdalabStatus.noiseSNR;
        save(TEMPFILE,'Ynoise', '-v7.3','SNR');
        
    end

end
