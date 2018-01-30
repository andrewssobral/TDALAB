function status=loadTensor()
    global loadedVars loadedFileName;
    global TDALABHOME;

    status=false;
    
    [fname, path] = uigetfile( '*.mat', 'Select a MAT-file',horzcat(TDALABHOME,filesep,'benchmark',filesep));
    if isequal(fname,0)
        status=false;
        return;
    end
    loadedFileName=horzcat(path,fname);
    
%     loadedVars=load(loadedFileName,'-mat');
    
    fvars=whos('-file',loadedFileName);

    if isempty(fvars)
        errordlg('The file does not contain any data.','Open error','modal');
        status=false;
    else
        loadedVars={fvars(:).name};  
        clsVars={fvars(:).class};
        bytes=[fvars(:).bytes];
        TDALABinitialization;
        h=findobj('tag','txtSelVar');
        set(h,'enable','on');
        varnames=loadedVars;
        h=findobj('tag','ppSelVar');
        set(h,'string',varnames,'enable','on','value',1);
        status=true;
        totalVarnames=numel(varnames);
        flag=true;
        for k=1:totalVarnames
            clsvar=clsVars{k};
            if strcmp(clsvar,'double')||strcmp(clsvar,'tensor')||strcmp(clsvar,'ttensor')||strcmp(clsvar,'ktensor')
                set(h,'value',k);
                flag=false;
                break;
            end
        end
        if flag
            h=errordlg('No acceptable variables found. Please select another file.','Load error','modal');
            waitfor(h);
            loadTensor;
        else
            h=allchild(0);
            oo = findall(h, 'Tag', 'TDALAB' );
            set(oo,'name',[fname,' :: ',defstr('TDALAB')]);
            selectVars;
        end
    end
end

