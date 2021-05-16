clc;
clear;
global KTYPE KSCALE size_training  C C_star data_flag size_up size_start C_end s
s = 0;
C=1;
C_star=1;
C_end=128;
KSCALE=0.5;  
temp = [];
temp1=0;
size_up=500;
for data=4
    data_flag=data;
    if data==1                  %w6a  sample:17188 unlabe:17000
        size_start=500;
    end
    if data==3                  %CodRNA  sample:59535 unlabel:59035
        size_start=500;
    end
    if data==2                  %text    sample:1946 unlabel:1800
        size_start=500;
    end
    if data==6                  %ijcnn1  sample:49990 unlabel:49790
        size_start=500;
    end
    if data==5                  %madelon  sample:2600 unlabel:2400
        size_start=500;
    end
    if data==4                  %usps     sample:2007 unlabel:1800
        size_start=3000;
    end
    if data==7                  %A9a  sample:32561 unlabel:32200
        size_start=500;
    end
    if data==8                  %Mushrooms  sample:8124 unlabel:7900
        size_start=500;
    end
    
    for k=1  
        if k==1
            KTYPE=6;
        end
        if k==2
            KTYPE=2;
        end
        if k==3
            KTYPE=1;
        end

        for j=1:1
            size_training=size_start+size_up*(j-1);
            for i=1:1
                [time,error,original_x,original_y,extended_x,extended_y,initial_solution,extended_real_y]=main(data_flag);
            end
        end
    end
end
% save('2_10000_5');