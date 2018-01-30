function [ h ] = plot3Dtensor(hax,Y,info)
%PLOT3DTENSOR Summary of this function goes here
%   Detailed explanation goes here

if ~isfield(info,'colormap') info.colormap='Jet'; end;
if ~isfield(info,'viewpoint') info.viewpoint=[45 25]; end;

Y=tennormalize(Y,inf);




plotfunc=info.plotfunc;

typeIplot={'contourf','contour','contour3'};
typeIIplot={'mesh','waterfall','stem3','surf','surfl'};

if strcmpi(class(Y),'ktensor')    
    subs=(1:length(Y.lambda))';
    subs=subs(:,ones(1,3));
    G=sptensor(subs,Y.lambda);
    Y=ttensor(tensor(G),Y.U);
end
Gsz=size(Y.core);
Ysz=size(Y);
if numel(Gsz)~=3
    error('Works for 3-way tensor ONLY.');
end

if isempty(hax)
    figure;
    hax=gca;
end
subplot(hax);
cla;

if strcmpi(info.core,'cube')
    cubeCore;
else
    [x y z]=meshgrid(0:(Gsz(2)+1)/Gsz(2):Gsz(2),0:1+1/Gsz(1):Gsz(1),0:1+1/Gsz(3):Gsz(3));
    slice(x,y,z,double(Y.core),[0 max(x(:))],[0,max(y(:))],[0,max(z(:))]);
    shading interp;
end
% colormap(info.colormap);

if strcmpi(class(Y),'double')
    return;
end

%% Continue to draw mode matrices
hold on;

L=1.5*max(Gsz)./Ysz;
space=1;
plotfc=str2func(plotfunc);

%% mode - 2
%% ####
c=2;
x0=Gsz(1)+space;
switch plotfunc
    case typeIplot
        [x y]=meshgrid(x0+[1:Ysz(2)]'.*L(c),0:(Gsz(2)+1)/Gsz(2):Gsz(2));
        plotfc(x,y,Y.U{2}');
    case typeIIplot
        [x y]=meshgrid(x0+[1:Ysz(2)]'.*L(c),0:(Gsz(2)+1)/Gsz(2):Gsz(2));
        plotfc(x,y,Y.U{2}');
    otherwise
        x=repmat([x0+[0:Ysz(2)-1]'.*L(c)],1,Gsz(2));
        y=bsxfun(@plus,Y.U{2}./2,[.5:Gsz(2)-.5]);
        z=zeros(Ysz(2),Gsz(2));
        plot3(x,y,z);
        bx=[x0,x0+(Ysz(2)-1).*L(c) x0+(Ysz(2)-1).*L(c) x0 x0]';
        by=[0 0 Gsz(2) Gsz(2) 0]';
        bz=zeros(5,1);
        plot3(bx,by,bz,'r-');
end
text((x0+Ysz(2).*L(c))/1.5,Gsz(2)*1.2,0,'{\bfA}_{2}');


%% ###########################  mode -1
c=1;
y0=Gsz(3)+space;
t1=hgtransform('parent',hax);
switch plotfunc
    case typeIplot
        [x y]=meshgrid(0:(Gsz(1)+1)/Gsz(1):Gsz(1),y0+[0:Ysz(1)-1]'.*L(c));
        [temp h1]=plotfc(x,y,Y.U{1});
        set(h1,'parent',t1);
    case typeIIplot
        [x y]=meshgrid(0:(Gsz(1)+1)/Gsz(1):Gsz(1),y0+[0:Ysz(1)-1]'.*L(c));
        h1=plotfc(x,y,Y.U{1});
        set(h1,'parent',t1);
    otherwise
        y=repmat(y0+[0:Ysz(1)-1]'.*L(c),1,Gsz(1));
        x=bsxfun(@plus,Y.U{1}./2,[.5:Gsz(1)-.5]);
        z=zeros(Ysz(1),Gsz(1));
        h1=plot3(x,y,z);
        set(h1,'parent',t1);
        by=[y0,y0+(Ysz(1)-1).*L(c) y0+(Ysz(1)-1).*L(c) y0 y0]';
        bx=[0 0 Gsz(1) Gsz(1) 0]';
        bz=zeros(5,1);
        hb1=plot3(bx,by,bz,'r-');set(hb1,'parent',t1);
end
ht1=text(Gsz(1)/3,y0+Ysz(1).*L(c)*1.1,0,'{\bfA}_{1}');
set(ht1,'parent',t1);

rot=makehgtform('translate',[0 Gsz(2) 0],'xrotate',pi/2);
set(t1,'matrix',rot);
% END of mode-1


%% ###########################  mode - 3
c=3;
y0=2*space;
t3=hgtransform('parent',hax);
switch plotfunc
    case typeIplot
        [x y]=meshgrid(0:(Gsz(3)+1)/Gsz(3):Gsz(3),-y0-[0:Ysz(3)-1]'.*L(c));
        [temp h3]=plotfc(x,y,Y.U{3});
        set(h3,'parent',t3);
    case typeIIplot
        [x y]=meshgrid(0:(Gsz(3)+1)/Gsz(3):Gsz(3),-y0-[0:Ysz(3)-1]'.*L(c));
        h3=plotfc(x,y,Y.U{3});
        set(h3,'parent',t3);
    otherwise
        y=repmat(-y0-[0:Ysz(3)-1]'.*L(c),1,Gsz(3));
        x=bsxfun(@plus,-Y.U{3}./2,[.5:Gsz(3)-.5]);
        z=zeros(Ysz(3),Gsz(3));
        h3=plot3(x,y,z);
        set(h3,'parent',t3);
        by=-[y0,y0+(Ysz(3)-1).*L(c) y0+(Ysz(3)-1).*L(c) y0 y0]';
        bx=[0 0 Gsz(3) Gsz(3) 0]';
        bz=zeros(5,1);
        hb3=plot3(bx,by,bz,'r-');set(hb3,'parent',t3);
end
ht3=text(Gsz(3)/3,-y0-Ysz(3).*L(c)*1.1-2,0,'{\bfA}_{3}');
set(ht3,'parent',t3);

rot=makehgtform('translate',[0 0 Gsz(3)],'yrotate',pi/2);
set(t3,'matrix',rot);
%   END of mode-3


text(Gsz(1)*1.1,Gsz(2)*1.1,Gsz(3)*0.75,'$\underline{\bf{G}}$','interp','latex');

title(info.title);
axis off;

view(info.viewpoint);
colormap(info.colormap);

%% ------------------------------------------------------------------------
    function cubeCore
        thres=0.01;
        G3=abs(double(Y.core));
        G3=G3/max(G3(:));
        I=size(G3);
        G3(G3<thres)=0;
        Xcube=0.5:I(1)-0.5;
        Ycube=0.5:I(2)-0.5;
        Zcube=0.5:I(3)-0.5;

        MX=max(Xcube);mX=min(Xcube);
        MY=max(Ycube);mY=min(Ycube);
        MZ=max(Zcube);mZ=min(Zcube);
        FV=cube(0.5*[(MX+mX),(MY+mY),(MZ+mZ)],0.5*[MX-mX,MY-mY,MZ-mZ]+0.5);
        % FV.facecolor='b';
        patch(FV,'Edgecolor','r','facealpha',0.1,'linewidth',1);

        for i=1:I(1)
        %     set(hwaitbarPatch,'XData',[0 100*i/I(1) 100*i/I(1) 0]);
            for j=1:I(2)
                for k=1:I(3)
                    if G3(i,j,k)>thres
                        FV=cube([Xcube(i),Ycube(j),Zcube(k)],G3(i,j,k)/2);
                        hc=patch(FV,'Edgecolor',[0.1 0.1 0.1]);
                    end
                end
            end
        end
    end
axis equal;
end

