function hf=showPIs(LEFTPOS,LEFTBOTTOM)
% A is a N-by-M matrix whose rows consists of the SIR values concerning mode-n
global SIRs NumOfMode NumOfComp tdalabStatus;
global XSpread YSpread ScreenWidth ScreenHeight  lFontSize;

if ~tdalabStatus.decomposed
    warndlg('Please run algorithm first.','Warning','modal');
    return;
end

if isempty(SIRs)
    errordlg('No SIRs information available.','Error','modal');
    return;
end

h=findobj('tag','tdvfigssir');
if ~isempty(h)
    set(h,'visible','on');
    return;
end
    

A=SIRs';

vp=3;
Width=120*XSpread;Height=40*YSpread;
allW=1.5*Width;
if nargin==0
    LEFTPOS=(ScreenWidth-allW)/2;
    LEFTBOTTOM=(ScreenHeight-Height)/2;
elseif nargin==1
    LEFTBOTTOM=(ScreenHeight-Height)/2;
end
hf=figure('Units','Characters','Resize','on','toolbar','figure',...
    'OuterPosition',[LEFTPOS,LEFTBOTTOM,Width,Height],...
    'name','Show SIRs - 3D bar','numbertitle','off','tag','tdvfigssir');
bkcolor=get(hf,'color');
pos=get(hf,'position');
Height=pos(4);Width=pos(3);
haxesBar=axes('parent',hf,'Units','Characters','position',[10*XSpread 4*YSpread Width-35*XSpread Height-6*YSpread],...
    'color',bkcolor,'fontsize',lFontSize);

hbar=bar3(A);
% 
%%%%%%%%%%%%%%%%%%%%%% visualization
Z = magic(max(NumOfComp));
for i = 1:length(hbar)
    zdata = ones(6*length(hbar),4);
    k = 1;
    for j = 0:6:(6*length(hbar)-6)
        zdata(j+1:j+6,:) = Z(k,i);
        k = k+1;
    end
    set(hbar(i),'Cdata',zdata)
end

shading interp
for i = 1:length(hbar)
    zdata = get(hbar(i),'Zdata');
    set(hbar(i),'Cdata',zdata)
    set(hbar,'EdgeColor','k')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end


title('SIRs of each component in each mode');
xlim([0.5 size(A,2)+0.5]);
ylim([0.5 size(A,1)+0.5]);
xlabel('Component');
ylabel('Mode');
zlabel('SIR (dB)');
zH=zlim;
zlim([0 zH(2)]);
view(vp);
colormap Cool
hcolorbar=colorbar('fontsize',lFontSize);



bkcolor=get(hf,'Color');
fkcolor='k';



shift=23*XSpread;
htvp=uicontrol('parent',hf,'Style','text','Units','Characters','Position',[Width-22*XSpread Height-3.3*YSpread 18*XSpread 2*YSpread],'string',...
    'ViewPoint:','backgroundcolor',bkcolor,'ForegroundColor',fkcolor,...
    'Fontsize',lFontSize,'Fontweight','bold');
heditvp=uicontrol('parent',hf,'Style','Edit','Units','Characters','Position',[Width-20*XSpread Height-6*YSpread 15*XSpread 2*YSpread],'value',...
    3,'backgroundcolor','white','ForegroundColor',fkcolor,'string','3',...
    'Fontsize',lFontSize,'Fontweight','bold',...
    'Tag','editvp');
hpbChange=uicontrol('parent',hf,'Style','pushbutton','Units','Characters','Position',[Width-20*XSpread Height-9*YSpread 15*XSpread 2*YSpread],'String',...
    'Change Now!','Tag','pbChange','callback',@changeView);
htbRotate=uicontrol('parent',hf,'Style','togglebutton','Units','Characters','Position',[Width-20*XSpread Height-12*YSpread 15*XSpread 2*YSpread],'String',...
    'Rotate3D','Tag','ptRotate','callback','rotate3d');
hpbHideCtrl=uicontrol('parent',hf,'Style','togglebutton','Units','Characters','Position',[Width-20*XSpread Height-15*YSpread 15*XSpread 2*YSpread],'String',...
    'Hide All','Tag','pbPIHideAll','callback',@hideAll);
hpbPICloseAll=uicontrol('parent',hf,'Style','togglebutton','Units','Characters','Position',[Width-20*XSpread Height-18*YSpread 15*XSpread 2*YSpread],'String',...
    'Close all','Tag','pbPICloseAll','callback',@closeAll);


pos=get(hf,'Outerposition');
Height=pos(4);

Width=pos(3);
space=YSpread/YSpread;
hspace=XSpread/XSpread;
figure('Units','Characters','Resize','on','toolbar','none',...
    'OuterPosition',[LEFTPOS+Width+hspace,LEFTBOTTOM+Height/2,Width/2,Height/2],...
    'name','Show SIRs-imagesc','numbertitle','off','tag','tdvfigssir');
imagesc(SIRs');
set(gca,'YTick',1:NumOfMode,'XTick',1:max(NumOfComp));
ylabel('Num. of modes');
xlabel('Num. of Comp.');
colorbar;



pmode={'-r.';'-bo';'-gp';'-ms';'-cv';'-yx'};
figure('Units','Characters','Resize','on','toolbar','none',...
    'OuterPosition',[LEFTPOS+hspace+Width,LEFTBOTTOM,Width/2,Height/2],...
    'name','Plot of SIRs','numbertitle','off','tag','tdvfigssir');
for n=1:NumOfMode
    plot(SIRs(:,n),pmode{n});
    hold on;
    axis tight;
    xlabel('mSIR(dB)');
    ylabel('Num. of modes.');
    grid('on');
    ls{n}=['Mode ',num2str(n)];
end
legend(ls);



    function changeView(hObject, eventdata, handles)
        heditvp=findobj('tag','editvp');
        z=str2num(get(heditvp,'String'));
        Lz=min(length(z),3);
        vp=z(1:Lz);
        if length(vp)>Lz
            vp(Lz+1:end)=[];
        end
        if (length(vp)==1)&&(vp(1)~=2)&&(vp(1)~=3)
            vp=3;
            msgbox('Error viewpoint values are specified! Please use 2 or 3 or x,y,z.','Input Error','error');
        end
        
        view(vp);
    end

    function hideAll(hObject,eventdata,handles)
%         ctrls=allchild(hf);
%         ctrls=setdiff(ctrls,[hbar hcolorbar get(hbar(1),'parent')]);
        set([htvp heditvp hpbChange htbRotate hpbHideCtrl hpbPICloseAll],'visible','off');
        set(hf,'units','normalized','toolbar','figure');
        set(get(hbar(1),'parent'),'units','normalized','position',[.15 .05 .7 .94]);
        set(hcolorbar,'units','normalized','position',[.9 .05 .05 .85]);
    end

    function closeAll(hObject, eventdata, handles)
        h=findobj('name','Plot of SIRs','-or',...
            'name','Show SIRs - 3D bar','-or',...
            'name','Show SIRs-imagesc');
        delete(h);
    end
end

