function saveWorkSpace
global Ycap Ynoise tdalabStatus algs algOptions  TDALABHOME loadedFileName;
global MCTimes MCSIRs SIRs elapsedTime MCelapsedTime fit MCfit TEMPFILE;

if ~tdalabStatus.decomposed
    errordlg('Decompose the tensor first before saving results.');
    return;
end
 

    h=allchild(0);
    oo = findall(h, 'Tag', 'TDALAB' );
    benchmark=get(oo,'name');
    benchmark=regexp(benchmark,' :: ','split');
    benchmark=strrep(benchmark{1},'.mat','');

    if tdalabStatus.advEvaluation
        [file,path] = uiputfile('*.mat','Save Workspace As',horzcat(TDALABHOME,filesep,'userdata',filesep,'MCExp_',benchmark,'_',date,'.mat'));
    else
        [file,path] = uiputfile('*.mat','Save Workspace As',horzcat(TDALABHOME,filesep,'userdata',filesep,'Exp_',benchmark,'_',date,'.mat'));
    end
    
    if isequal(file,0)||isequal(path,0)
        errordlg('The workspace is not saved.','Failed to save');
        return;
    end

if isempty(file)
    errordlg('Cannot save the file.','File saving error.');
    return;
end


    dataFile=loadedFileName;
    Yname=findobj(allchild(oo),'tag','ppSelVar');
    Y=get(Yname,'string');
    Y=Y(get(Yname,'value'));
    Y=Y{1};
    

    varList=cell(1,7);
    if tdalabStatus.advEvaluation
        algNames={algs(tdalabStatus.algIndex).name};
        allparas={algOptions{tdalabStatus.algIndex}};
        for idx=1:numel(algNames)
            algParas.(algNames{idx})=allparas{idx};
        end
        varlist={'tdalabStatus','algNames','algParas','MCTimes','MCelapsedTime','MCSIRs','MCfit','Y','dataFile'};
    else
        algParas=algOptions{tdalabStatus.algIndex};
        algNames=algs(tdalabStatus.algIndex).name;
        varlist={'tdalabStatus','algNames','algParas','elapsedTime','SIRs','fit','Y','Ycap','dataFile'};
    end

    if ~isinf(tdalabStatus.noiseSNR)
        if exist(TEMPFILE,'file')~=0
            load(TEMPFILE,'Ynoise','SNR');
            noiseSNR=SNR;
            varlist{end+1}='Ynoise';
            varlist{end+1}='noiseSNR';
        end
    end

    timestamp=horzcat('Saved at ',datestr(now));
    varlist{end+1}='timestamp';

    %% construct values
    for n=1:numel(varlist)
        if ~isempty(varlist{n})
            Exp.(varlist{n})=eval(varlist{n});
        end
    end
    

    save(horzcat(path,file),'Exp','-v7.3');



end