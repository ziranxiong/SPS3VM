function [x,y]=read_data(flag)
global unlabel_size 
    if flag==1
        load w6a;
        unlabel_size=17000;
    end
    if flag==3
        load codrna;
        unlabel_size=59035;
%         unlabel_size=9035;
    end
    if flag==2
        load letter;
        unlabel_size=19950;
    end
    if flag==6
        load ijcnn1;
        unlabel_size=49790;
    end
    if flag==5
        load dota2;
        unlabel_size=92600;
    end
    if flag==4
        load usps;
        %unlabeled_index=find(y<=4);
        %labeled_index=[1:length(y)];
        %labeled_index(unlabeled_index)=[];
        %y(unlabeled_index) = -1;
        %y(labeled_index) = 1;
        %save uspst x y;
        unlabel_size=7241;
    end
    if flag==7
        load a9a;
        unlabel_size=32200;
    end
    if flag==8
        load mushrooms;
        unlabel_size=7900;
    end
    x=zscore(x);
end
