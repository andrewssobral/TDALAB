function waveplot
global NumOfComp  NumOfMode  sirs  Ycap Y ;
global   tdalabStatus;
% global hwaitbarPatch hwaitbar;
global    lFontSize  defaultFontName  ;
algIndex=tdalabStatus.algIndex;
Category=tdalabStatus.model;
if isempty(Ycap)
    Ycap=Y;
%     warndlg('Please run algorithm first.','Warning','modal');
%     return;
end

NumOfCompcap=cell2mat(cellfun(@(x) size(x,2),Ycap.U,'uni',false));

Acap=Ycap.U;
Gcap=Ycap.lambda;

if strcmpi(tdalabStatus.inputType,'tensor')
    A=[];G=[];PARAFAC_lambda=[];
else
    A=Y.U;
    if strcmpi(tdalabStatus.model,'CP')
        PARAFAC_lambda=Y.lambda;
    else
        G=Y.core;
    end
end

r=0.42;
figUp=0.84;
figAllSpace=0.00;
hfigAllInOne=[];
WSpace=0.32;
subFigsWidth=WSpace-figAllSpace;
hG3d=[];hGest3d=[];hSsig=[];hSest=[];
hfigAllInOne=[];

NumOfPlots=10;CurrentMode=1;
oo=findobj('name','Wave Plot');
if ~isempty(oo)
    delete(oo);
end

if (~strcmpi(tdalabStatus.inputType,'tensor'))&&strcmp(Category,'CP')
    % correct the signs
    %% Permute the results
    for n=1:NumOfMode
        [sirs(:,n) orders(:,n)]=CalcSIR(A{n},Acap{n});
    end
    [temp, minPos]=max(min(sirs));
    
    for n=1:NumOfMode
        A{n}=A{n}(:,orders(:,minPos));
    end
    pos=zeros(1,max(NumOfComp));
    flag=pos;
    for n=1:NumOfMode
        pos=(sign(diag((A{n})'*Acap{n}))<0)&(sirs(:,n)>10);
        if n~=NumOfMode
            flag(pos)=1-flag(pos);
            Acap{n}(:,pos)=-Acap{n}(:,pos);
        else
            Acap{n}(:,flag>0)=-Acap{n}(:,flag>0);
        end
    end
end

% if strcmpi(tdalabStatus.inputType,'tensor')
%     mainW=0.5;
% else
    mainW=0.97;
% end

hmainWave=figure('Units','normalized','outerposition',[0.5-mainW/2 0.89 mainW 0.1],'toolbar','none',...
    'menubar','none','tag','mainWave','name','Wave plot','numbertitle','off');
set(hmainWave,'CloseRequestFcn',@closeWave);


%% components
bkcolor=get(hmainWave,'color');
hgpMode=uibuttongroup('parent',hmainWave,'Units','normalized','position',[0.001 0.5 0.998 0.49],...
    'backgroundcolor',bkcolor,'bordertype','line','visible','on','SelectionChangeFcn',@selMode,...
    'tag','gpWaveSelMode');
set(hgpMode,'borderwidth',1,'highlightcolor','k');
step=0.999/NumOfMode;
for n=1:NumOfMode
    htbModen=uicontrol('parent',hgpMode,'Units','normalized','position',[0.001+(n-1)*step 0.01 step 0.98],...
        'string',['Mode' num2str(n)],'style','togglebutton','fontsize',lFontSize,'fontname',defaultFontName,...
        'userdata',num2str(n),'backgroundcolor',[0.52 0.53 0.67]);
end
hgpNum=uibuttongroup('parent',hmainWave,'Units','normalized','position',[0.001 0.001 0.998 0.49],...
    'backgroundcolor',bkcolor,'bordertype','line','visible','on','SelectionChangeFcn',@selNum);
set(hgpMode,'borderwidth',1,'highlightcolor','k');
step=0.8/11;
htxt=uicontrol('parent',hmainWave,'Units','normalized','position',[0.001 0.14 0.1 0.25],...
    'backgroundcolor',bkcolor,'visible','on','style','text','string','Num. of sigs per page:',...
    'fontsize',lFontSize,'fontname',defaultFontName);
for i=1:10
    h=uicontrol('parent',hgpNum,'Units','normalized','position',[0.101+(i-1)*step 0.01 step 0.98],...
        'string',num2str(i),'style','radiobutton','fontsize',lFontSize,'fontname',defaultFontName,'backgroundcolor',bkcolor);
end

h=uicontrol('parent',hmainWave,'Units','normalized','position',[0.94 0.001 0.059 0.49],'style','pushbutton',...
    'string','Close all','callback',@closeWave);
htpAllInOne=uicontrol('parent',hmainWave,'Units','normalized','position',[0.84 0.001 0.059 0.49],'style','togglebutton',...
    'string','All in one','callback',@drawAllInOne,'backgroundcolor',bkcolor,'tag','tbWavesAllInOne');

set(hgpNum,'SelectedObject',[]);
drawWaves(1);
set(htpAllInOne,'value',0);


    function selMode(source,eventdata)
        if strcmp(eventdata.EventName,'SelectionChanged')
            str=get(eventdata.NewValue,'userdata');
            CurrentMode=str2num(str);
            drawWaves(1);
            if get(findobj('tag','tbWavesAllInOne'),'value')==1
                drawAllInOne;
            end
            OldMode=str2num(get(eventdata.OldValue,'userdata'));
            h=findobj('name',['Source signals of mode ' num2str(OldMode)],...
                '-or','name',['Estimated signals of mode ' num2str(OldMode)]);
            delete(h);
        end
    end
    function selNum(source,eventdata)
        NumOfPlots=str2num(get(eventdata.NewValue,'string'));
        drawWaves(0);
    end


    function drawWaves(modeChanged)
        isDrawAll=get(htpAllInOne,'value');
        if strcmpi(tdalabStatus.inputType,'tensor')
            if isDrawAll==0
                subFigsWidth=(0.97-figAllSpace)/2;
            else
                subFigsWidth=(0.97-2*figAllSpace)/3;
            end
            if strcmp(Category,'CP')
                subs=(1:NumOfCompcap(1))';subs=subs(:,ones(1,NumOfMode));
                Gcap=sptensor(subs,Ycap.lambda);
            end
            param = struct( 'signal', (Acap{CurrentMode})', ...
                    'standard',0,...
                    'x_limit',[],...
                    'x_title',[],...
                    'num', NumOfPlots, ...
                    'stitle', ['Estimated signals of mode ' num2str(CurrentMode)], ...
                    'sname', 'y', ...  % changed from 'x'
                    'title_first_plot', [], ...
                    'sposition', [ 0.015 0.05+r subFigsWidth figUp-r ], ...  %[ 0.01 0.05 0.48 0.37 ]
                    'sunits', 'normalized', ...
                    'window_number', 0 );
                hSest=icalab_srplot( param );
                if modeChanged==1
                    G3d.name='Structure of tucker core G [Estimation]';
                    G3d.position=[0.015+subFigsWidth 0.05+r subFigsWidth figUp-r];
                    hGest3d=visual3D(Gcap,CurrentMode,G3d);
                end
        else % A is not empty
            if isDrawAll==0
                subFigsWidth=(0.97-figAllSpace)/2;
            else
                subFigsWidth=(0.97-2*figAllSpace)/3;
            end
            if strcmp(Category,'CP')
                subs=(1:NumOfComp(1))';subs=subs(:,ones(1,NumOfMode));
                G=sptensor(subs,PARAFAC_lambda);
                if strcmp(class(Ycap),'ttensor')
                    Gcap=Ycap.core;
                else
                    subs=(1:NumOfCompcap(1))';subs=subs(:,ones(1,NumOfMode));
                    Gcap=sptensor(subs,Ycap.lambda);
                end
            end
            param = struct( 'signal', (A{CurrentMode})', ...
                'standard',0,...
                'x_limit',[],...
                'x_title',[],...
                'num', NumOfPlots, ...
                'stitle', ['Source signals of mode ' num2str(CurrentMode)], ...
                'sname', 'y', ...  % changed from 'x'
                'title_first_plot', [], ...
                'sposition', [ 0.015 0.05+r subFigsWidth figUp-r ], ...  %[ 0.01 0.05 0.48 0.37 ]
                'sunits', 'normalized', ...
                'window_number', 0 );
            
            hSsig=icalab_srplot( param );
            param.signal=(Acap{CurrentMode})';param.stitle=horzcat('Estimated signals of mode ', num2str(CurrentMode));
            %                 param.sposition=[ 0.343 0.05+r 0.315 figUp-r ];
            param.sposition=[ 0.015 0.05 subFigsWidth r ];
            hSest=icalab_srplot( param );
            if modeChanged==1
                G3d.name='Structure of tucker core G [Original]';
                %                     G3d.position=[0.67 0.05+r 0.315 figUp-r];
                G3d.position=[0.015+subFigsWidth 0.05+r subFigsWidth figUp-r];
                hG3d=visual3D(G,CurrentMode,G3d);
                
                G3d.name='Structure of tucker core G [Estimation]';
                G3d.position=[0.015+subFigsWidth 0.05 subFigsWidth r];
                hGest3d=visual3D(Gcap,CurrentMode,G3d);
            end
        end
    end
drawnow

    function drawAllInOne(hObject,eventdata,handles)
        if get(findobj('tag','tbWavesAllInOne'),'value')==1
            if strcmpi(tdalabStatus.inputType,'tensor')
                if subFigsWidth>0.4
                    subFigsWidth=(0.97-2*figAllSpace)/3;
                    set(hSest,'outerposition',[0.015 0.05+r subFigsWidth figUp-r]);
                    set(hGest3d,'outerposition',[ 0.015+2*(subFigsWidth+figAllSpace) 0.05+r subFigsWidth figUp-r]);
                end
                hfigAllInOne=figure('Unit','Normalized','outerposition',[ 0.015+(subFigsWidth+figAllSpace) 0.05+r subFigsWidth figUp-r],'Name',horzcat('Estimated signals of mode ',num2str(CurrentMode)),...
                    'numbertitle','off');
                hold on;
                for n=1:size(Acap{CurrentMode},2)
                    plot(Acap{CurrentMode}(:,n),'-.');
                end
            else % A is not empty
                if subFigsWidth>0.4
                    subFigsWidth=(0.97-2*figAllSpace)/3;
                    set(hSsig,'outerposition',[0.015 0.05+r subFigsWidth figUp-r]);
                    set(hSest,'outerposition',[0.015 0.05 subFigsWidth r]);
                    set(hG3d,'outerposition',[0.015+2*(subFigsWidth+figAllSpace) 0.05+r subFigsWidth figUp-r]);
                    set(hGest3d,'outerposition',[ 0.015+2*(subFigsWidth+figAllSpace) 0.05 subFigsWidth r]);
                end
                hfigAllInOne(1)=figure('Unit','Normalized','outerposition',[0.015+subFigsWidth+figAllSpace 0.05 subFigsWidth r],'Name',horzcat('Estimated signals of mode ',num2str(CurrentMode)),...
                    'numbertitle','off');
                hold on;
                hfigAllInOne(2)=figure('Unit','Normalized','outerposition',[ 0.015+subFigsWidth+figAllSpace 0.05+r subFigsWidth figUp-r],'Name',horzcat('Source signals of mode ',num2str(CurrentMode)),...
                    'numbertitle','off');
                hold on;
                figure(hfigAllInOne(1));
                for n=1:NumOfCompcap(CurrentMode)
                    plot(Acap{CurrentMode}(:,n),'-.');
                end
                figure(hfigAllInOne(2));
                for n=1:NumOfComp(CurrentMode)
                    plot(A{CurrentMode}(:,n),'-.');
                end
            end
        else % delete allInOne
            if subFigsWidth<0.4
                subFigsWidth=(0.97-figAllSpace)/2;
                if ~strcmpi(tdalabStatus.inputType,'tensor')
                    set(hSsig,'outerposition',[0.015 0.05+r subFigsWidth figUp-r]);
                    set(hSest,'outerposition',[0.015 0.05 subFigsWidth r]);
                    set(hG3d,'outerposition',[0.015+(subFigsWidth+figAllSpace) 0.05+r subFigsWidth figUp-r]);
                    set(hGest3d,'outerposition',[ 0.015+(subFigsWidth+figAllSpace) 0.05 subFigsWidth r]);
                else
                    set(hSest,'outerposition',[0.015 0.05+r subFigsWidth figUp-r]);
                    set(hGest3d,'outerposition',[ 0.015+(subFigsWidth+figAllSpace) 0.05+r subFigsWidth figUp-r]);
                end
            end
            if ~isempty(hfigAllInOne)&&ishandle(hfigAllInOne)
                delete(hfigAllInOne);
            end
        end
    end


    function closeWave(hObject,eventdata,handles)
        h=([hSsig hSest hG3d hGest3d hfigAllInOne]);
        h=h(ishandle(h));
        if ~isempty(h)
            delete(h);
        end
        h=findobj('name','Wave plot');
        delete(h);
    end



end

