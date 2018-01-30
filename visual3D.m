function haxes=visual3D(modes,plotOpts,haxes)
%% plotOpts.{G,title,opos,tag}
%  modes={[1 3],[2],[4]}, as an example.
if ~isfield(plotOpts,'units')
    plotOpts.units='normalized';
end

if nargin==2
    haxes=[];
end

thres=0.1; % g is ignored if g/max(plotOpts.signal)< thres 
% global hwaitbarPatch hwaitbar;
MaxL=20;
Tol=1e-3;
NumOfComp=size(plotOpts.signal);NumOfMode=numel(NumOfComp);
if NumOfMode==2
    %% added in the future;
%     info.signal=double(plotOpts.signal);
%     info.plotfunc='mplot';
%     mplotmatrix(haxes,info);
    h=figure('Units',plotOpts.units,'outerposition',plotOpts.opos,'name',plotOpts.title,...
        'numbertitle','off','tag',plotOpts.tag);
    contour(double(plotOpts.signal));
%     colormap gray;
    axis tight;
    return;
end


if numel(unique(cell2mat(modes)))~=NumOfMode
    error('Incompatible partition of modes.');
end
xmode=modes{1};
ymode=modes{2};
zmode=modes{3};
xl=prod(NumOfComp(xmode));
yl=prod(NumOfComp(ymode));
zl=prod(NumOfComp(zmode));
G3=double(permute(plotOpts.signal,[xmode ymode zmode]));
G3=reshape(G3,xl,yl,zl);
% G3=permute(G3,[3 2 1]);
G3=abs(G3);

I=size(G3);
if I(1)>MaxL
    r=MaxL/I(1);
    G3=G3(floor(1:r:I(1)),:,:);
end
if I(2)>MaxL
    r=MaxL/I(2);
    G3=G3(:,floor(1:r:I(2)),:);
end
if I(3)>MaxL
    r=MaxL/I(3);    
    G3=G3(:,:,floor(1:r:I(3)));
end
I=size(G3);
R=max(G3(:))/2;
G3(G3<R*thres*2)=0;
X=2*R*(1:I(1));
Y=2*R*(1:I(2));
Z=2*R*(1:I(3));

if isempty(haxes)
    h=figure('Units',plotOpts.units,'outerposition',plotOpts.opos,'name',plotOpts.title,...
        'numbertitle','off','tag',plotOpts.tag);
    haxes=gca;
else
    set(get(haxes,'parent'),'units',plotOpts.units);
    cla;
end

MX=max(X);mX=min(X);
MY=max(Y);mY=min(Y);
MZ=max(Z);mZ=min(Z);
FV=cube(0.5*[(MX+mX),(MY+mY),(MZ+mZ)],0.5*[MX-mX,MY-mY,MZ-mZ]+R);
FV.facecolor='b';
patch(FV,'Edgecolor','k','facealpha',0.1,'linewidth',1);

for i=1:I(1)
%     fprintf('i=%d\n',i);
%     set(hwaitbarPatch,'XData',[0 100*i/I(1) 100*i/I(1) 0]);
    for j=1:I(2)
        for k=1:I(3)
            if G3(i,j,k)>Tol
                FV=cube([X(i),Y(j),Z(k)],G3(i,j,k)/2);
                patch(FV,'Edgecolor',[0.1 0.1 0.1],'linestyle','-');
            end
        end
    end
end
% close(hwaitbar);

xlabel(horzcat('Dim -',num2str(xmode,' [%d]')));
ylabel(horzcat('Dim -',num2str(ymode,' [%d]')));
zlabel(horzcat('Dim -',num2str(zmode,' [%d]')));
set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
title(plotOpts.title);
% axis equal;
axis image;
view(60,30);