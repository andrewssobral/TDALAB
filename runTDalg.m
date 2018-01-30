function runTDalg
global algOptions algs tdalabStatus paraList hmainTDALAB elapsedTime;
global Y Ycap NumOfComp NumOfMode SIRs fit;
TDModelstr=defstr('TDModel');
% Run the algorithm
tdalab('hide');
commandwindow;
pause(0.1);

if isempty(tdalabStatus.algIndex)
    errordlg('Please select an algorithm first.','Setting error','modal');
    return;
else
    if tdalabStatus.advEvaluation        
        choice=questdlg('Configuration for Monte-Carlo test is detected. Run Monte-Carlo test?','Confirm','Monte Carlo Run','Normal Run','Normal Run');
        if strcmpi(choice,'Monte Carlo Run')
            MCTDRun;
            return;
        else
            tdalabStatus.algIndex=tdalabStatus.algIndex(1);
            tdalabStatus.advEvaluation=false;
        end
    end
end

fprintf('[TDALAB] Initializing ... \n');
%% initialize
SIRs=[];Ycap=[];Ynoise=[];

algIndex=tdalabStatus.algIndex(1);
tdalg=str2func(algs(algIndex).name);

fprintf('[TDALAB] Checking the parameters ...\n');
%% options
if isempty(algOptions{algIndex})
    algOptions{algIndex}=struct;
end

%% specify the number of components automatically
if isfield(algOptions{algIndex},'NumOfComp')
    if strcmpi(tdalabStatus.model,'CP')
        algOptions{algIndex}.NumOfComp=algOptions{algIndex}.NumOfComp(1);
    elseif strcmpi(tdalabStatus.model,'Tucker')&&numel(algOptions{algIndex}.NumOfComp)==1
        algOptions{algIndex}.NumOfComp=repmat(algOptions{algIndex}.NumOfComp,1,NumOfMode);
    end
end

fnames=fieldnames(algOptions{algIndex});
for n=1:numel(fnames)
    if isempty(algOptions{algIndex}.(fnames{n}))
        algOptions{algIndex}=rmfield(algOptions{algIndex},fnames{n});
    end
end

%% load observations Ynoise [tdalabtemp.mat]
global TEMPFILE;
fprintf('[TDALAB] Preparing the data for decomposition ... \n');
if get(findall(allchild(0),'tag','ckToTensor'),'value')==true
    if exist(TEMPFILE,'file')~=0
        load(TEMPFILE,'Ynoise','SNR');
    else
        Ynoise=[];SNR=rand;
    end
    if isempty(Ynoise)||(SNR~=tdalabStatus.noiseSNR)  %% need Ynoise
        if isinf(tdalabStatus.noiseSNR)
                fprintf('[TDALAB] Generating the full tensor. This may cost a few minutes depending on the problem size. Please wait ...\n');
                Ynoise=tensor(Y);

            if strcmpi(class(Ynoise),'double')
                Ynoise=tensor(Ynoise);
            end
            SNR=tdalabStatus.noiseSNR;
            save(TEMPFILE,'Ynoise','-v7.3','SNR');
        else
            h=errordlg('[TDALAB] Noise data not found. Please check your noise settings.','Invalid input','modal');
            uiwait(h);
            tdalab;
            return;
        end
    end
else
    Ynoise=Y;
end


%% prepare for running
dispRunInformation;

switch tdalabStatus.model
    case {'CP','Tucker','PMF'}
        fprintf('[TDALAB] Algorithm is running. Please wait...\n');
        ts=tic;
        if strcmpi(tdalabStatus.model,'PMF')
            Ynoise=double(Ynoise);
            try
                [y w]=tdalg(Ynoise,algOptions{algIndex});
            catch ME
                toc(ts);
                algerr;
                return;
            end
            Ycap=ktensor(ones(size(y,2),1),{y,w'});
        else
            if numel(size(Ynoise))<3
                fprintf(2,'Number of modes is less than 3. Please use PMF/BSS model.\n');
                tdalab;
                return;
            end
            try
                Ycap=tdalg(Ynoise,algOptions{algIndex});
            catch ME
                toc(ts);
                algerr;
                return;
            end
        end
        elapsedTime=toc(ts);
        set(hmainTDALAB,'HandleVisibility','callback','visible','on');
        commandwindow;
        
        clsYcap=class(Ycap);
        if strcmpi(clsYcap,'ttensor')
            Ycap=ttnormalize(Ycap,0);
        elseif strcmpi(clsYcap,'ktensor')
            Ycap=ktnormalize(Ycap,0);
        end
        %%
        if ~isempty(Ycap)
            fprintf('## Calculating the fitting error. This may cost several minutes depending the data size ...\n');
            [fit res]=fitness(Y,Ycap);
            fprintf('## Elapsed time: %f sesconds.  Fit: %f\n',elapsedTime,fit);
            tdalabStatus.decomposed=true;
            if (~strcmpi(tdalabStatus.inputType,'tensor'))
                %         SIRs = cellfun(@CalcSIR,Ycap.U(:),Y.U(:),'uni',0)';
                SIRs=NaN(max(NumOfComp),NumOfMode);
                if numel(NumOfComp)==1
                    NumOfComp=NumOfComp(ones(1,NumOfMode));
                end
                for n=1:NumOfMode
                    if size(Y.U{n},2)==size(Ycap.U{n},2)
                        SIRs(1:NumOfComp(n),n)=CalcSIR(Y.U{n},Ycap.U{n});
                    end
                end
                if ~all(isnan(SIRs))
                    SIRs=real(SIRs);
                    CPdispSIRs(SIRs);
                end
            end 
            
            if strcmpi(tdalabStatus.model,'Tucker')    
                h=findobj('tag','gbExtraConstraints');
                hsel=get(h,'SelectedObject');
                consType=get(hsel,'string');
                if ~isempty(consType)&&~strcmpi(consType,'none')
                    fprintf('[TDALAB] Imposing constraints on the mode matrices ... \n');                
                    t=refineTucker;
                    elapsedTime=elapsedTime+t;
                    fprintf('[TDALAB] Total time consumption: %f. \n',elapsedTime);
                end
            end
        else
            tdalabStatus.decomposed=false;
            fprintf('[TDALAB] Mission fails. Please carefully check your settings.\n');
        end
    case {'Parafac2'}
        fprintf('[TDALAB] Algorithm is running. Please wait...\n');
        ts=tic;
        try
            Ycap=tdalg(Ynoise,algOptions{algIndex});
        catch ME
            toc(ts);
            algerr;
            return;
        end
        elapsedTime=toc(ts);
        set(hmainTDALAB,'HandleVisibility','callback','visible','on');
        commandwindow;
        fprintf('[TDALAB] Mission complete. Elapsed elapsedTime is %f seconds.\n',elapsedTime);
        if ~isempty(Ycap)
            tdalabStatus.decomposed=true;
            if ~strcmpi(tdalabStatus.inputType,'tensor')
                SIRs(1:NumOfComp(1),1)=CalcSIR(Y.U{1},Ycap.A);
                SIRs(1:NumOfComp(1),3)=CalcSIR(Y.U{3},Ycap.C);
                SIRs(1:NumOfComp(1),2)=CalcSIR(Y.U{2},Ycap.P{1}*Ycap.H);
                CPdispSIRs(SIRs);
            end
        else
            tdalabStatus.decomposed=false;
            fprintf('[TDALAB] Mission fails. Please carefully check your settings.\n');
        end
    case 'BCD'
        Ynoise=tensor(Ynoise);
        fprintf('[TDALAB] Algorithm is running. Please wait...\n');
            ts=tic;
        try
            Ycap=tdalg(Ynoise,algOptions{algIndex});
        catch ME
            toc(ts);
            algerr;
            return;
        end
        elapsedTime=toc(ts);
        set(hmainTDALAB,'HandleVisibility','callback','visible','on');
        commandwindow;
        if ~isempty(Ycap)
            tdalabStatus.decomposed=true;
            fprintf('Mission complete. Elapsed elapsedTime is %f seconds.\n',elapsedTime);
            disp('--------------------------------------------------------------');
            disp('Struct: Ycap.');
            disp(Ycap);
        else
            tdalabStatus.decomposed=false;
        end
    otherwise
        fprintf(2,'[TDALAB] Unsuported tensor decomposition model.\n');
        return;
end

if tdalabStatus.decomposed
    set(allchild(findobj('tag','pnlOutputAnalysis')),'enable','on');
else
    set(allchild(findobj('tag','pnlOutputAnalysis')),'enable','off');
end

commandwindow;
fprintf('\n========================== END ==========================\n\n');
fprintf('\n[TDALAB] Type or click <a href="matlab: tdalab">tdalab</a> to return.\n\n');


    function algerr
        tdalab;
        commandwindow;   
        fprintf(2,'[TDALAB] Error occured during running the [%s] algorithm:\n',dispname(algs(algIndex).details));
        fprintf(2,horzcat('[TDALAB] ',ME.message));      
        for id=numel(ME.stack):-1:1
            disp(ME.stack(id));
        end
%         disp('[TDALAB] Type or click <a href = "matlab: tdalab">tdalab</a> to return.');
        fprintf('\n======================= <a href="matlab: tdalab">TDALAB</a> =======================\n\n');   
    end

end
