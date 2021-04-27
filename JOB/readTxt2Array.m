function [array] =  readTxt2Array(filepath,n)
    file = fopen(filepath,'rt');
    if(file == -1)
        disp(['Error opening the file:',filepath]);
    end
    array = cell(1,n);
    count =1;
    while 1
        nextline = fgetl(file);
        if ~ischar(nextline)
            break;
        end
        a=sscanf(nextline,'%s%s');
   
        %class(a);
        array(count)={a};
        count=count+1;
    end
end