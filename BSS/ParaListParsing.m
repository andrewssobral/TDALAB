function paraList=ParaListParsing(paraList,paraTypeList,indices)

%% work for matlab2010a and above
% for i=1:numel(indices)
%     currAlg=indices(i);
%     names=fieldnames(paraTypeList{currAlg});
%     for c=1:numel(names)
%         currna=names{c};
%         switch paraTypeList{currAlg}.(currna)
%             case paraType.paraDouble
%                 if ~isfloat(paraList{currAlg}.(currna))
%                     paraList{currAlg}.(currna)=str2num(paraList{currAlg}.(currna));
%                 end
%             case paraType.paraMenu
%                 p=regexp(paraList{currAlg}.(currna),'\|','split');
%                 paraList{currAlg}.(currna)=p{1};
%             case paraType.paraTF
%                 paraList{currAlg}.(currna)=paraList{currAlg}.(currna);
%             otherwise
%                 error('Unsupported parameter type. -- in [PMFsetOpts].');
%         end
%     end
% end

for i=1:numel(indices)
    currAlg=indices(i);
    names=fieldnames(paraTypeList{currAlg});
    for c=1:numel(names)
        currna=names{c};
        if paraTypeList{currAlg}.(currna)==paraType.paraDouble
            if ~isfloat(paraList{currAlg}.(currna))
                paraList{currAlg}.(currna)=str2num(paraList{currAlg}.(currna));
            end
        elseif paraTypeList{currAlg}.(currna)==paraType.paraString
            paraList{currAlg}.(currna)=char(paraList{currAlg}.(currna));
        elseif paraTypeList{currAlg}.(currna)==paraType.paraMenu
            p=regexp(paraList{currAlg}.(currna),'\|','split');
            paraList{currAlg}.(currna)=p{1};
        elseif paraTypeList{currAlg}.(currna)==paraType.paraTF
            paraList{currAlg}.(currna)=logical(paraList{currAlg}.(currna));
        else
            error('Unsupported parameter type. -- in [PMFsetOpts].');
        end
    end
end
end