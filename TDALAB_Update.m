function TDALAB_Update
updateURL='http://bsp.brain.riken.jp/TDALAB';
% verURL='file:///C:/zhougx/matlabcodes/TDALAB_Ver0.6xx/tdalab_ver.ini';
verURL='http://bsp.brain.riken.jp/~zhougx/tdalab/tdalab_ver.ini';

    c_ver=tdalab_version;
    getver=urlread(verURL);
    getver=getver(1:100);
    getver=regexp(getver,'\s*\n','split');
    getver=regexp(getver{2},' ','split');
    getver=cellfun(@str2num,getver,'uni',false);
    flag=cellfun(@isempty,getver);
    getver=getver(flag==0);
    if isempty(getver)
        errordlg('No version information available. Please visit this function later.','Error');
        return;
    else
        getver=getver{1};
    end
    if isequal(getver,c_ver)
        helpdlg('Your version is the latest. Thank you for checking.','Update');
    else
        choice=questdlg(horzcat('New version TDALAB Ver. ',num2str(getver),' is available. Do you want to download it now?'),'Update now?','Yes','No','Yes');
        switch choice
            case 'Yes'
                web(updateURL,'-browser');
            case 'No'
                return;
        end
    end
    
end