function h=imagematrix(info)
% by Guoxu ZHOU
% imagesc/contourf/surfc
definfo=struct('signal',randn(2,100),'title','Image','imagefunc','imagesc','colormap','hsv','width',[],'cols',[],'stitle','s',...
    'tag','imagematrix','units','normalized','opos',[],'viewpoint',[60 30]);
if ~exist('info','var')
    info=struct;
end
info.signal=real(info.signal);
sz=size(info.signal);
info=scanparam(definfo,info);
MODE=numel(sz);
if MODE==3
    M=sz(3);
    
    if M>20
        display('Too many to plot. Only the first 20 slices will be displayed.');
        M=20;
    end
    
    info.width=sz(2);
    height=sz(1);
elseif numel(sz)==2
    M=sz(1);T=sz(2);
    
   if M>20
        display('Too many to plot. Only the first 20 slices will be displayed.');
        M=20;
    end
    
    height=floor(T/info.width);
    info.signal=info.signal(:,1:height*info.width);
end
if isempty(info.cols)
    if M<=5
        info.cols=M;
    else
        info.cols=ceil(sqrt(M));
    end
%     K=sqrt(M);
%     for k=floor(K):-1:3
%         if ~rem(M,k)
%             info.cols=k;
%             break;
%         end
%     end
%     if isempty(info.cols)
%         info.cols=ceil(K);
%     end
end

h=findobj('tag',info.tag,'-and','name',info.title);
if isempty(h);
    h=figure('name',info.title,'units',info.units,'numbertitle','off','outerposition',info.opos,'tag',info.tag);
else
    figure(h);
    clf;
    if isempty(info.opos)
        move(h,'center');
    else
        set(h,'units',info.units,'outerposition',info.opos);
    end
end
set(gca,'fontname','Times New Roman');
cla;



space=0.01;
rows=ceil(M/info.cols);
axw=(1-(info.cols-1)*space)/info.cols;
axh=(.85-(rows)*space)/rows;
cla;
for c=1:info.cols
    for r=1:rows
        k=(r-1)*info.cols+c;
        if k>M
            break;
        end
        hsubplot=subplot('position',[(c-1)*(axw+space) .85-r*axh-r*space axw axh]);
        if MODE==3
            sk=info.signal(:,:,k);
        else
            sk=reshpe(info.signal(k,:),info.width,height);
        end
        switch info.imagefunc
            case 'imagesc'
                imagesc(sk);
            case 'contourf'
                contourf(flipud(sk));
            case 'surf'
                surf(sk,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
                view(info.viewpoint);
            otherwise
                fc=str2func(info.imagefunc);
                fc(sk);
        end
%         if ~isempty(info.stitle)
%             th=title(horzcat(info.stitle,'-',num2str(k)));
%             set(th,'FontSize',8);
%         end
        if r==1
            th=title(horzcat('Col.',num2str(c)));
            set(th,'FontUnit','point','Fontsize',8);
        end
        axis off;
    end % row
end % cols
colormap(info.colormap);
end
