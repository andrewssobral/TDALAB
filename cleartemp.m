function cleartemp()
global Ycap sirs tdalabStatus TEMPFILE;
Ynoise=[];SNR=[];
if ~isempty(TEMPFILE)
    save(TEMPFILE,'Ynoise','SNR');
end
Ycap=[];
sirs=[];
tdalabStatus.decomposed=false;