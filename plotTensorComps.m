function [ h ] = plotTensorComps( Y,opts )

defopts=struct('plotfunc','plot','position',[]);
if ~exist('opts','var')
    opts=struct();
end
[plotfunc,fig_pos]=scanparam(defopts,opts);

Ycls=class(Y);
tensorName=inputname(1);

switch Ycls
    case 'ktensor'
        NumOfComp=length(Y.lambda);
        NumOfMode=numel(size(Y));
        NumOfComp=NumOfComp(ones(1,NumOfMode));
        h=showComps;
    case 'ttensor'
        NumOfComp=size(Y.core);
        NumOfMode=numel(NumOfComp);
        h=showComps;
    otherwise
        h=[];
        error('Unsupported data type.');
end

    function h=showComps
        tops=0.1;bots=0.1;Left=0.1;Rs=0.05;
        
        h=figure('Unit','inches','Position',[0. 0. 3.5 2],'name',tensorName,'NumberTitle','off');
        set(h,'Unit','Normalized');
        movegui(h,'center');
        
        
        height=(1-tops-bots)/(NumOfMode);
          
        
        Top=1-tops;
        for n=1:NumOfMode
            bottom=Top-n*height;
            width=(1-Left-Rs)/NumOfComp(n);
            for c=1:NumOfComp(n)
                subplot('position',[Left+(c-1)*width,bottom,width,height]);
                drawComp(Y.U{n}(:,c),plotfunc);
                set(gca,'XTickLabel','','YTickLabel','');
                if c==1
                    ylabel([num2str(n)]);
                end
            end
        end 
    end
    function drawComp(s,plotfunc)  % draw a signal by using plotfunc
        switch plotfunc
            case 'plot'
                plot(s);
                axis tight;
                smin=min(s);smax=max(s);
                ylim([smin-abs(smin)*0.1 smax+abs(smax)*0.1]);
                
            otherwise 
                error('Unsupported plot function.');
        end        
    end


end

