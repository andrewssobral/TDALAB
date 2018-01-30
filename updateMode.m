function updateMode
global NumOfComp NumOfMode Y tdalabStatus;
if isempty(tdalabStatus.inputType)
    errordlg('No tensor/ktensor/ttensor found.','No input error','modal');
    return;
elseif strcmpi(tdalabStatus.inputType,'tensor')
    return;
end

h=findobj('tag','ppMode');
n=get(h,'value');

[s selchn]=selectchn((Y.U{n})');
hrbTucker=findobj('tag','rbTucker');

ind=1:size(Y.U{n},2);
selchn=ind(selchn==1);

if ~isempty(s)
    if strcmp(class(s),'double')
        s=s';
        if strcmp(tdalabStatus.model,'CP')
            if size(s,2)~=NumOfComp(1)
                errordlg('All the mode matrices must have the same columns in CP model.','Valid selection','modal');
                return;
            else
                Y.U{n}=s;
                cleartemp;
            end
        else % Tucker model
            A=Y.U;G=double(Y.core);A{n}=s;
            sz=size(Y);
            ind=cell(1,NumOfMode);
            for i=1:NumOfMode
                ind{i}=1:NumOfComp(i);
            end
            ind{n}=selchn;
            G=G(ind{:});
            Y=ttensor(tensor(G),A);
            cleartemp;
            NumOfComp(n)=size(Y.U{n},2);
            tdalabStatus.inputType='ttensor';
        end
        cleartemp;
    elseif strcmp(s,'delete')&&NumOfMode>3
        if strcmp(tdalabStatus.inputType,'ktensor')
            lambda=Y.lambda;A=Y.U;
            lambda=lambda(selchn);A=A{selchn};
            Y=ktensor(lambda,A);
            cleartemp;
        elseif strcmp(tdalabStatus.inputType,'ttensor')
            G=double(Y.core);G=sum(G,n);G=squeeze(G);A=Y.U;
            A(n)=[];
            Y=ttensor(tensor(G),A);
            cleartemp;
        else
            errordlg('Unknown error during updating the tensor.','Unknown error','modal');
        end
        NumOfMode=NumOfMode-1;
        NumOfComp(n)=[];
    else
        errordlg(['Fail to update mode-',num2str(n)],'Update error','modal');
        return;
    end
    updateUI;

end
