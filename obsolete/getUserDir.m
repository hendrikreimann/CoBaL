



function user_dir = getUserDir
    if ispc
        error('Function getUserDir not implemented yet for Windows systems. Fix this!')
    else
        user_dir = getenv('HOME');
    end
end