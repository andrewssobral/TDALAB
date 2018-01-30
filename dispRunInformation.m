function dispRunInformation
global tdalabStatus NumOfMode NumOfComp;
global algs algOptions Y;
algIndex=tdalabStatus.algIndex;

title=['\n\n===============   Tensor decomposition [',tdalabStatus.model,' model]   =============== \n'];
fprintf(title);

    fprintf(horzcat('* Input type           : ',tdalabStatus.inputType,'\n'));
        str=horzcat('* Dimension of tensor : %d',repmat(' x %d',1,NumOfMode-1),'\n');
    fprintf(str,size(Y));
        str=horzcat('* Number of Components : ',num2str(NumOfComp,' %d'),'.');
    disp(str);
    if isinf(tdalabStatus.noiseSNR)
        fprintf('* Noise free.\n');
    else
        fprintf('* Additive noise       : [%s] SNR=%ddB, Sparsity: %d%%.\n',tdalabStatus.noiseType,tdalabStatus.noiseSNR,tdalabStatus.noiseSparsity);
    end
    fprintf(horzcat('* Algorithm            : ',dispname(algs(algIndex).details),'\n'));
    fprintf('* Options:\n');
    disp(algOptions{algIndex});
