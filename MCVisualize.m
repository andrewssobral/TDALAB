function MCVisualize

global algs tdalabStatus NumOfMode;
global MCfit MCSIRs MCelapsedTime;

NAlgs=numel(tdalabStatus.algIndex);
algNameList=cellfun(@(x) regexprep(x,'_','\\_'),{algs(tdalabStatus.algIndex).name},'uni',false);
%% fitness
fig_MCfit=figure('units','inches','position',[0 0 3.5 2.3],'name','Fit','numbertitle','off','tag','tdvfigsMCfit','visible','off','RendererMode','manual','render','painters');
movegui(fig_MCfit,'center');
MCpos=get(fig_MCfit,'position');
MCpos(1)=MCpos(1)-1.75;MCpos(2)=MCpos(2)+2.3;
set(fig_MCfit,'position',MCpos);
outMCpos=get(fig_MCfit,'outerposition');
set(fig_MCfit,'units','normalized','visible','on');

lines=plot(MCfit','*-');

mks=defstr('Marker');
nmks=numel(mks);
%         lines=get(gca,'child');
cmaps=cool(numel(lines));
for idx=1:numel(lines)
    set(lines(idx),'Marker',mks(rem(idx-1,nmks)+1),'Color',cmaps(idx,:));
end

grid on;
legend((strrep(algNameList,'call\_','')),'Location','southeast');
ylim([min(min(MCfit(:)),-0.1) 1.1]);
xlabel('Run');
ylabel('Fitting error');
axis tight;
ylim([0 1.05]);

fig_sirs=[];

%% sir
if (~strcmpi(tdalabStatus.inputType,'tensor'))
    fig_sirs=figure('units','inches','outerposition',[outMCpos(1)+outMCpos(3) outMCpos(2) outMCpos(3) outMCpos(4)],'name','mean SIRs','numbertitle','off','tag','tdvfigsMCfit','RendererMode','manual','render','painters');
    set(gcf,'units','normalized');
    % x - mode y - algs z - sirs
    hbar=bar3(mean(MCSIRs,3));
    hbar0=get(hbar(1),'parent');
    for n=1:NumOfMode
        modestr{n}=['mode - ',num2str(n)];
    end
    set(gca,'XTick',1:NumOfMode,'XTickLabel',modestr,'YTick',1:NAlgs,'YTickLabel',(strrep({algs(tdalabStatus.algIndex).name},'call_','')));
    zlabel('meanSIRs (dB)');
    %%%%%%%%%%%%%%%%%%%%%% visualization
    Z = magic(max(NumOfMode,NAlgs));
    for i = 1:length(hbar)
        zdata = ones(6*length(hbar),4);
        k = 1;
        for j = 0:6:(6*length(hbar)-6)
            zdata(j+1:j+6,:) = Z(k,i);
            k = k+1;
        end
        set(hbar(i),'Cdata',zdata)
    end
    
    shading interp
    for i = 1:length(hbar)
        zdata = get(hbar(i),'Zdata');
        set(hbar(i),'Cdata',zdata)
        set(hbar,'EdgeColor','k')
    end
    hcb=colorbar;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end
    set(hbar0,'position',[0.13 .08 .62 0.9],'FontSize',8);
    set(hcb,'position',[0.89 .1 .04 0.8],'FontSize',8);
    axis tight;
end


%% time
fig_time=figure('units','inches','outerposition',[outMCpos(1)+outMCpos(3)/2 outMCpos(2)-outMCpos(4) outMCpos(3) outMCpos(4)],'name','Elapsed time','numbertitle','off','tag','tdvfigsMCfit','RendererMode','manual','render','painters','visible','off');
bar(mean(MCelapsedTime,2));
ylabel('Elapsed time (s)');
axis tight;
pos=axis(gca);
axis([pos(1)-.5 pos(2)+.5 pos(3) pos(4)+.2]);
set(gca,'XTick',1:NAlgs,'XTickLabel',(strrep({algs(tdalabStatus.algIndex).name},'call_','')),'FontSize',8);
if isempty(fig_sirs)
    set(fig_time,'outerposition',[outMCpos(1)+outMCpos(3) outMCpos(2) outMCpos(3) outMCpos(4)]);
end
set(fig_time,'visible','on');


txtMCResults;


    function txtMCResults
        nWidth=20;
        i=1;
        fill=' ';
        algNames=fill(ones(NAlgs,nWidth));
        for algIndex=tdalabStatus.algIndex
            name=algs(algIndex).name;
            lname=length(name);
            if lname>nWidth
                algNames(i,1:nWidth-3)=name(1:17);
                algNames(i,nWidth-2:nWidth)='...';
            else
                algNames(i,1:lname)=name;
            end
            i=i+1;
        end
        fprintf('\n------------------------------- Results ------------------------------\n');
        title=' Algorithms        Time(s)   Fit      ';
        for n=1:NumOfMode
            title=[title,'  meanSIR(',num2str(n),')'];
        end
        fprintf([title,'\n']);
        formatstr=['%-18s%6.1f    %6.2f',repmat('    %8.4f',1,NumOfMode),'\n'];
        mtimes=mean(MCelapsedTime,2);
        mMCfits=mean(MCfit,2);
        mMCSIRs=mean(MCSIRs,3);
        for i=1:NAlgs
            fprintf(formatstr,algNames(i,:),mtimes(i),mMCfits(i),mMCSIRs(i,:));
            i=i+1;
        end
        for algIndex=tdalabStatus.algIndex
            d=dispname(algs(algIndex).details);
            fprintf('* [%-15s] -->   [%s]\n',algs(algIndex).name,d);
        end
    end

end %% MCVisualize