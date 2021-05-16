function [per]=getAcc(y1,y2)
    count = 0;
    for i=1:length(y1)
        if y1(i) == y2(i)
            count = count + 1;
        end
    end
    per = count/length(y1);
end
