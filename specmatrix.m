function specmatrix(haxes,info)
% by Guoxu ZHOU

definfo=struct('signal',sin(0:.1:2*pi),'window','hamming','winlength',128,'noverlap',8,'tag','specmatrix','opos',[0.01 0.01 0.98 0.98],'units','normalized',...
    'title','Spectrogram analysis');
if ~exist('info','var')
    info=struct;
end
info=scanparam(definfo,info);

if isempty(haxes);
    figure('name',info.title,'units',info.units,'numbertitle','off','outerposition',info.opos,'tag',info.tag);
    haxes=gca;
else
    set(get(haxes,'parents'),'units',info.units);
end
axes(haxes);
set(gca,'fontname','Times New Roman');
cla;

[M T]=size(info.signal);
cols=floor(M*5/7);
rows=ceil(M/cols);
space=0.1;
axw=(1-(cols+0.5)*space)/cols;
axh=(1-(rows+0.5)*space)/rows;

info.winlength=2^min(info.winlength,floor(log2(T/4)));
if info.noverlap>info.winlength
    info.noverlap=info.noverlap/2;
end
winfunc=str2func(info.window);
win=window(winfunc,info.winlength);
for c=1:cols
    for r=1:rows
        k=(r-1)*cols+c;
        if k<=M
            subplot('position',[(c-1)*(axw+space)+space 1-r*axh-(r-0.5)*space axw axh]);
            spectrogram(info.signal(k,:),win,info.noverlap,info.winlength);
        end
    end
end       

end