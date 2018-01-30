function str=dispname(str)
str=regexp(str,'</.*>','split');
str=regexp(str{1},'<.*>','split');
str(cellfun(@isempty,str))=[];
str=char(str);
end