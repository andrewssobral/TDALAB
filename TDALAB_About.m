function TDALAB_About
global ScreenWidth ScreenHeight XSpread YSpread lFontSize defaultFontName;
Width=85*XSpread;
Height=25*YSpread;
hmainAbout=figure('Units','Characters','position',[0.5*ScreenWidth-0.5*Width,ScreenHeight/2, Width, Height],...
    'menu','none','toolbar','none','resize','on','name','About...','numbertitle','off','color','white',...
    'tag','mainAbout','windowstyle','modal');
set(hmainAbout,'Units','pixels');
%pos=get(hmainAbout,'position');
set(hmainAbout,'Units','Characters');
bkcolor=get(hmainAbout,'color');
load('smallbsi.mat');
hlog=axes('parent',hmainAbout,'Units','normalized','position',[0.0 0.1 0.4*YSpread/XSpread 0.4],'color',bkcolor);
image(logofig);
axis image;axis off;

set(hlog,'Units','normalized');
pos=get(hlog,'position');
left=pos(1)+pos(3)+0.01;
width=1-left-.1;

verstr=horzcat('TDALAB Ver.',num2str(tdalab_version,'%3.1f'));
uicontrol('parent',hmainAbout,'Units','normalized','position',[0.01, 0.83, 0.98 0.15],'fontweight','normal',...
    'style','text','string',verstr,'fontunits','normalized','fontsize',.6,'Max',8,'Min',1,'fontname',defaultFontName,'backgroundcolor','white');


str={'Guoxu ZHOU and Andrzej CICHOCKI';'';...
    'Laboratory for Advanced Brain Signal Processing';'';...
    'Brain Science Institute, RIKEN, Japan'};
uicontrol('parent',hmainAbout,'Units','normalized','position',[0.01, 0.52, 0.98 0.3],'fontweight','normal',...
    'style','text','string',str,'fontunits','normalized','fontsize',.15,'Max',8,'Min',1,'fontname',defaultFontName,'backgroundcolor','white');

uicontrol('parent',hmainAbout,'Units','normalized','position',[left, 0.38, width 0.1],...
    'style','pushbutton','string','http://bsp.brain.riken.jp','fontsize',lFontSize,'foregroundcolor','b','enable','inactive',...
    'fontname',defaultFontName,'buttondownfcn','web(''http://bsp.brain.riken.jp'','' -browser'')',...
    'TooltipString','Visit the homepage.');

uicontrol('parent',hmainAbout,'Units','normalized','position',[left, 0.24, width 0.1],...
    'style','pushbutton','string','Feedback & bug report','fontsize',lFontSize,'foregroundcolor','b','enable','inactive',...
    'fontname',defaultFontName,'buttondownfcn','web(''mailto:zhouguoxu@brain.riken.jp'')');

uicontrol('parent',hmainAbout,'Units','normalized','position',[left+width/2-0.1 0.04 0.2 0.12],...
    'style','pushbutton','string','Close','callback','h=findobj(''tag'',''mainAbout'');delete(h);','fontsize',lFontSize);
end