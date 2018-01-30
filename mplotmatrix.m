function mplotmatrix(haxes,info)
% by Guoxu ZHOU
definfo=struct('signal',rand(2,100),'xlabel','Index','ylabel','\itt','zlabel','\it{s}({\itt})','viewpoint',[25 45],'title','plotmatrix','linewidth',2,...
    'linestyle','','plotfunc','plot3','colormap','hsv','tag','plotmatrix','opos',[0.01 0.01 0.98 0.98],'units','normalized','legend','none');
if ~exist('info','var')
    info=struct;
end
info=scanparam(definfo,info);

if isempty(haxes);
    figure('name',info.title,'units',info.units,'numbertitle','off','outerposition',info.opos,'tag',info.tag);
    haxes=gca;
else
    set(get(haxes,'parent'),'units',info.units);
end
axes(haxes);
set(gca,'fontname','Times New Roman');
cla;

info.signal=real(info.signal);
[M T]=size(info.signal);

for m=1:M
%     complabel{m}=['Comp. ' num2str(m)];
    complabel{m}=[num2str(m)];
end

if min(M,T)>50
    imagesc(info.signal);
    axis equal;
    axis tight;
    colormap gray;
    title(info.title);
    return;
end

switch info.plotfunc
    case 'plot'
        view(2);
        plot(info.signal',info.linestyle,'linewidth',info.linewidth');
        grid on;
        xlabel(info.ylabel);
        ylabel(info.zlabel);
        axis tight;
        if ~strcmpi(info.legend,'none')
            c=regexp(info.legend,'\s*\|\s*','split');
            c=strtrim(c);
            c=c(~cellfun('isempty', c));
            legend(c);
        end
    case 'mplot'
        view(2);
        plot(bsxfun(@plus,datanormalize(info.signal',inf)./2,1:M),info.linestyle,'linewidth',info.linewidth');
        grid on;
        xlabel(info.ylabel);
        ylabel('Components');
        set(gca,'YTick',1:M);
        axis tight;
        
    case 'plot3'
        view(3);
        for m=1:M
            plot3(m(ones(1,T)),1:T,info.signal(m,:),info.linestyle,'linewidth',info.linewidth);
            hold on;
        end
        axis tight;
        L=axis;        
        fv.FaceColor='interp';
        fv.Faces=[1 2 3 4];
        fv.EdgeColor='none';
        for m=1:M
            fv.Vertices=[m,L(3),L(5);m,L(4),L(5);m,L(4),L(6);m,L(3),L(6)];
            fv.facevertexcdata=rand(length(fv.Vertices),3);
            fv.FaceAlpha=max(0.85-m/4,0.2);
            patch(fv);
        end
        hold off;
%         xlabel(info.xlabel);
        ylabel(info.ylabel);
        zlabel(info.zlabel);
        set(gca,'XTick',1:M,'XTickLabel',complabel);
        xlabel('Components');
        grid on;
        view(info.viewpoint);
    case 'stem3'
        view(3);
        stem3(info.signal,'o-');
        axis tight;
        xlabel(info.ylabel);
%         ylabel(info.ylabel);
        zlabel(info.zlabel);
        grid on;
        view(info.viewpoint);
        set(gca,'YTick',1:M,'YTickLabel',complabel);
    case 'ribbon'
        view(3);
        hrib=ribbon(info.signal');
        set(hrib,'LineStyle','none');
%         get(hrib(1))
        axis tight;
%         xlabel(info.xlabel);
        ylabel(info.ylabel);
        zlabel(info.zlabel);
        grid on;
        set(gca,'XTick',1:M,'XTickLabel',complabel);
        view(info.viewpoint);
    case 'waterfall'
        h=waterfall(info.signal);
        axis tight;
        set(h,'facevertexcdata',rand(length(get(h,'Vertices')),3),'FaceAlpha',0.4,'FaceColor','flat','linewidth',info.linewidth,...
            'Edgecolor','none');
        colormap(info.colormap);
        hold on;
        for m=1:M
            plot3(1:T,m(ones(1,T)),info.signal(m,:),'linewidth',info.linewidth);
            hold on;
        end
        axis tight;
%         ylabel(info.xlabel);
        xlabel(info.ylabel);
        zlabel(info.zlabel);
        grid on;
        view(info.viewpoint);
        set(gca,'YTick',1:M,'YTickLabel',complabel);
        ylabel('Components');
end
title(info.title);
end


