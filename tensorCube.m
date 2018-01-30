function h=tensorCube(hax,Y,opts)
defopts=struct('colormap','Jet','dimorder',[1 2 3],'position',[],'xlabel','1','ylabel','2','zlabel','index',...
    'direction','y');
if ~exist('opts','var')
    opts=struct();
end
[colormap,dimorder,pos,xstr,ystr,zstr,mode]=scanparam(defopts,opts);


NumOfMode=numel(size(Y));
order=1:NumOfMode;
order(dimorder)=[];
order=[dimorder order];
Y=permute(Y,order);
Ydim=size(Y);
Y=reshape(Y,[Ydim(1) Ydim(2) prod(Ydim(3:end))]);


if isempty(hax)
    h=figure;
    hax=gca;
else
    subplot(hax);
    h=gcf;
end
cla;

cubeCoreY;
axis tight;
view(3);

    function cubeCoreZ        
        G3=abs(double(Y));
        G3=G3/max(G3(:));
        I=size(G3);
                
        thres=0.01;
        G3(G3<thres)=0;
        Xcube=0.5:I(1)-0.5;
        Ycube=0.5:I(2)-0.5;
        Zcube=0.5:I(3)-0.5;

        MX=max(Xcube);mX=min(Xcube);
        MY=max(Ycube);mY=min(Ycube);
        MZ=max(Zcube);mZ=min(Zcube);
        
        %% z indices...        
        NZ=prod(Ydim(4:end));
        rzi=[(MZ-mZ)+1]/NZ;
        for n=1:NZ
            FV= cube(0.5*[(MX+mX),(MY+mY),2*(n-.5)*rzi],0.5*[MX-mX,MY-mY,rzi-1]+0.5); 
            patch(FV,'Edgecolor','g','facealpha',0.1,'linewidth',1,'linestyle',':');
        end
        
        FV=cube(0.5*[(MX+mX),(MY+mY),(MZ+mZ)],0.5*[MX-mX,MY-mY,MZ-mZ]+0.5);
        % FV.facecolor='b';
        patch(FV,'Edgecolor','r','facealpha',0.1,'linewidth',1.5);
        xlabel(xstr);ylabel(ystr);

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
    set(gca,'XTickLabel','','YTickLabel','','ZTick',Ydim(3):Ydim(3):NZ*Ydim(3),'ZTickLabel',1:NZ);
    end

    function cubeCoreY
        Y=permute(Y,[1 3 2]);
        G3=abs(double(Y));
        G3=G3/max(G3(:));
        I=size(G3);
                
        thres=0.01;
        G3(G3<thres)=0;
        Xcube=0.5:I(1)-0.5;
        Ycube=0.5:I(2)-0.5;
        Zcube=0.5:I(3)-0.5;

        MX=max(Xcube);mX=min(Xcube);
        MY=max(Ycube);mY=min(Ycube);
        MZ=max(Zcube);mZ=min(Zcube);
        
        %% z indices...        
        NZ=prod(Ydim(4:end));
        ryi=(MY-mY+1)/NZ;
        for n=1:NZ
            FV= cube(0.5*[(MX+mX),2*(n-.5)*ryi,(MZ+mZ)],0.5*[MX-mX,ryi-1,MZ-mZ]+0.5); 
            if rem(n,2)
                patch(FV,'Edgecolor','g','facealpha',0.3,'linewidth',1,'linestyle',':');
            else                
                patch(FV,'Edgecolor','g','facealpha',0.1,'linewidth',1,'linestyle',':');
            end
        end
        
        FV=cube(0.5*[(MX+mX),(MY+mY),(MZ+mZ)],0.5*[MX-mX,MY-mY,MZ-mZ]+0.5);
        % FV.facecolor='b';
        patch(FV,'Edgecolor','r','facealpha',0.1,'linewidth',1.5);
        xlabel(xstr);zlabel(ystr);ylabel(zstr);

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
    set(gca,'XTickLabel','','ZTickLabel','','YTick',Ydim(3):Ydim(3):NZ*Ydim(3),'YTickLabel',1:NZ);
    end

end