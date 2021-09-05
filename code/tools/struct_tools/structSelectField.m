function strOut = structSelectField(strIn,selectV)
% Apply a selection on all variables in a struct which match a particular     

    fnames = fieldnames(strIn);
   
    strOut = rmfield(strIn,setdiff(fnames,selectV));

end