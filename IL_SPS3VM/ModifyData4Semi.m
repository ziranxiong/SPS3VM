function [y]=ModifyData4Semi(y)
    global unlabel_size 
    length1=length(y);
    res = [];
    while length(res) ~= unlabel_size
        %tmp_res = unique(randint(unlabel_size,1,[1,length1]));
        tmp_res = unique(randi([1,length1],unlabel_size,1));
        tmp_res = sort(tmp_res);
        res = [res;tmp_res];
        res = unique(res);
        if length(res) >= unlabel_size
            res = res(1:unlabel_size);
            res = sort(res);
        end
    end
    y(res)=0;
end