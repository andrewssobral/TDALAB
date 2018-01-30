function TDALABconfig


global TDALABPATH;
TDALABPATH=genpath(pwd);
addpath(TDALABPATH);

global TDALABHOME TEMPFILE;
p=which('tdalab');
cd(p(1:end-8));
TDALABHOME=pwd;
TEMPFILE=horzcat(TDALABHOME,filesep,'userdata',filesep,'tdalabtemp.mat');

global XSpread YSpread ScreenWidth ScreenHeight Spread;
oldUnits=get(0,'Units');

set(0,'Units','characters');
spt=get(0,'MonitorPositions');
spt=spt(1,:);

set(0,'Units','characters');
WinSize=get(0,'ScreenSize');
XSpread=1;
YSpread=1;
Spread=min(XSpread,YSpread);

global PixelPt;
PixelPt=spt(3)/WinSize(3);

set(0,'Units','characters');
WinSize=get(0,'ScreenSize');
ScreenWidth=WinSize(3);
ScreenHeight=WinSize(4);
set(0,'Units',oldUnits);

%% UIFontSize -- normalized
global  nFontSize defaultFontName sFontSize lFontSize varFontName;
nFontSize=.8;
sFontSize=8; % fixed small font size
lFontSize=11; % fixed large font size
defaultFontName='Times New Roman';
varFontName=defaultFontName;

%% UIControlSize  -- characters
global defaultCtrlHeight defaultCtrlWidth defaultVSpace defaultHSpace;
defaultCtrlHeight=1.5*YSpread;
defaultCtrlWidth=30*XSpread;
defaultVSpace=YSpread;
defaultHSpace=XSpread;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% colormap
global backGroundColor colorMap pbBackGroundColor;
pbBackGroundColor=[0.9412 0.9412 0.9412];
backGroundColor = [0.831372549019608 0.815686274509804 0.784313725490196];                
colorMap = [ 1    0    0
                    0    1    0
                    0    0    1
                    1    1    0
                    0    1    1
                    1    0    1 ];


                

                 