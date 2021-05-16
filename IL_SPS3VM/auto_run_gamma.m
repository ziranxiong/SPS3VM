clc;
clear;
global KTYPE KSCALE size_training  C C_star data_flag size_up size_start C_end s
s = 0;
C=128;
C_star=0.1;
C_end=C;
KSCALE=0;  
sum = 0;
sum1 = 0; 
temp = [];
temp1=0;
size_up=10000;
for data=3:3
    data_flag=data;
    if data==2                  %letter    sample:1946 unlabel:1800
        size_start=500;
    end
    if data==3                  %CodRNA  sample:59535 unlabel:59035
        size_start=1000;
    end
    if data==4                  %usps     sample:2007 unlabel:1800
        size_start=500;
    end
    if data==5                  %dota2  sample:2600 unlabel:2400
        size_start=500;
    end
    if data==6                  %ijcnn1  sample:49990 unlabel:49790
        size_start=500;
    end
    if data==7                  %A9a  sample:32561 unlabel:32200
        size_start=500;
    end
    if data==8                  %Mushrooms  sample:8124 unlabel:7900
        size_start=1500;
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
        error_list = [];
        for j=1:1
            size_training=size_start;
            for i=1:10
                KSCALE = KSCALE + 0.2;
                [time,error,original_x,original_y,extended_x,extended_y,initial_solution,extended_real_y]=select_gamma(data_flag);
                error_list = [error_list; error];
                fprintf('i = %d\n', i);
            end
        end
    end
    save([num2str(data), '_Gamma'], 'error_list')
end
