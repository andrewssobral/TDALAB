function [ p ] = tdalab_version()
    ver=1.1;
    if nargout==0
        fprintf('TDALAB Version %3.1f\n',ver);
    else
        p=ver;
    end
end

