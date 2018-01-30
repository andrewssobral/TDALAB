function CPdispSIRs(sirs)
%%function CPdispSIRs(sirs)
% sirs is a matrix
% column j denotes the averaged sirs of each component in mode-j
global NumOfMode NumOfComp tdalabStatus;

    fprintf('## SIRs of each Comp.:\n');
    LNOC=max(NumOfComp);
    hstr='          ';
    for n=1:NumOfMode
        hstr=[hstr,'    Mode ',num2str(n)];
    end
    hstr=[hstr,'\n'];
    fprintf(hstr);
    k1str=horzcat('Comp. %d   ',repmat('%10.4f',1,LNOC));
    k10str=horzcat('Comp. %d  ',repmat('%10.4f',1,LNOC));
    for k=1:max(NumOfComp)
        if k<10
            fprintf(k1str,k,sirs(k,:));
            fprintf('\n');
        else
            fprintf(k10str,k,sirs(k,:));
            fprintf('\n');
        end
    end
    if strcmp(tdalabStatus.model,'CP')
        ms=mean(sirs);
    else
        for n=1:NumOfMode
            msn=sirs(:,n);msn(isnan(msn))=[];
            ms(n)=mean(msn);
        end
    end
    fprintf(['Avg.      ',repmat('%10.4f',1,LNOC),'\n'],ms);
    if strcmp(tdalabStatus.model,'Tucker')
        fprintf('\n*''NaN'': means that this component does not exist.\n');
    end
    fprintf('\n');
end