function [ h ] = visualizeClusters( Feas,est_labels,gnd_labels)
%VISUALIZECLUSTERS Summary of this function goes here
%   Detailed explanation goes here
if nargin==2
    gnd_labels=[];
end
nfea=size(Feas,2);
nfea;
h=[];
if nfea~=2&&nfea~=3
    error('The number of columns of Feas must be 2 or 3.');
    return;
end
if isempty(gnd_labels)
    cls=unique(est_labels);
else
    cls=unique(gnd_labels);
end
ncls=numel(cls);

h=figure('units','inch','position',[0 0 3.5 3.5],'visible','off','tag','tdvfigs_cls');

colors=lines(numel(est_labels));

clsscatter(est_labels,'.');
if ~isempty(gnd_labels)
    clsscatter(gnd_labels,'o');
end
axis tight;
if nfea==3
    view(3);
end
movegui(h,'center');
set(h,'visible','on');

    function clsscatter(labels,mkr)
        for n=1:ncls
            flag=(labels==n);
            if nfea==2
                plot(Feas(flag,1),Feas(flag,2),mkr,'color',colors(n,:));
                hold on;
            else
                plot3(Feas(flag,1),Feas(flag,2),Feas(flag,3),mkr,'color',colors(n,:));
                hold on;
            end
        end
    end

end
