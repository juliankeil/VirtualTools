function[]= vt_save_mat(mat,filename)
% based on save_mat by Stephan Moratti
% usage: save_mat(mat,filename);
% mat =  datamatrix
% filename = filename (full path)
% This function will generate a Tab-delimited file with data rounded to 4 decimals
%
% Julian Keil 2016
%%
[lines,columns]=size(mat);

fid = fopen(filename,'w');

if iscell(mat) == 1
    
    for l=1:lines
        for c = 1:columns
            if isnumeric(mat{l,c})
                fprintf(fid, '%.4f', mat{l,c});
            else
                fprintf(fid,mat{l,c});
            end
            if c<columns % don't print tab for rightmost column
                fprintf(fid,'\t');
            end
        end
        fprintf(fid,'\n');
    end
    
else    
    for l=1:lines
        for c = 1:columns
            fprintf(fid,'%.4f ',mat(l,c));
            if c<columns % don't print tab for rightmost column
                fprintf(fid,'\t');
            end
        end
        fprintf(fid,'\n');
    end
end
fclose(fid);
