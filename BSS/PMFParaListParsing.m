function PMFActParas=PMFParaListParsing(indices,PMFActParas,paraList,paraTypeList)
NumOfAlgs=numel(indices);
if nargin==1
    [paraList algs paraTypeList]=PMFalgInit();
    PMFActParas=cell(NumOfAlgs,1);
    for i=1:NumOfAlgs
        PMFActParas{i}=struct();
    end
end
for i=1:NumOfAlgs
    currAlg=indices(i);
    names=fieldnames(paraTypeList{currAlg});
    for c=1:numel(names)
        currna=names{c};
        if paraTypeList{currAlg}.(currna)==paraType.paraDouble
                if ~isfield(PMFActParas{i},currna)||~isfloat(PMFActParas{i}.(currna))
                    if isfloat(paraList{currAlg}.(currna))
                        PMFActParas{i}.(currna)=paraList{currAlg}.(currna);
                    else
                        PMFActParas{i}.(currna)=str2num(paraList{currAlg}.(currna));
                    end
                end
        elseif paraTypeList{currAlg}.(currna)==paraType.paraString
                p=regexp(paraList{currAlg}.(currna),'\|','split');
                if ~isfield(PMFActParas{i},currna)
                    PMFActParas{i}.(currna)=p{1};
                else
                    if ~find(strcmpi(PMFActParas{i}.(currna),p))
                        PMFActParas{i}.(currna)=p{1};
                    end
                end
        elseif paraTypeList{currAlg}.(currna)==paraType.paraTF
            if ~isfield(PMFActParas{i},currna)
                PMFActParas{i}.(currna)=paraList{currAlg}.(currna);
            end
        else
                error('Unsupported parameter type. -- in [PMFsetOpts].');
        end
    end    
end
end