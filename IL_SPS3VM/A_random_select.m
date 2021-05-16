function [training_x, training_y, test_x, test_y]=A_random_select(x, y)
% ���ѡȡ����
% �����ݼ�����Ϊ ���Լ�(test) & ѵ����(training)
% ���Լ���test_x, test_y
% ѵ������training_x, training_y
% ʹ�� size_training ����ѵ�����Ĵ�С
global size_training;
local_flag=1;
while(local_flag)
    training_x=[];
    training_y=[];
    local_length=length(x(:,1));
    test_x=x;
    test_y=y;
    len1=0;
    while len1<size_training
        % �������һ�� 1 * size_training-len1 �ľ��󣬾����е�������[1, length-len1]֮��
        local_rand= randi([1,local_length-len1],1,size_training-len1);
        % ȥ�������е��ظ�����
        local_rand=unique(local_rand);
        len1=len1+length(local_rand);
        local_x=test_x(local_rand,:);
        local_y=test_y(local_rand);
        training_x=[training_x; local_x];
        training_y=[training_y; local_y];
        test_x(local_rand, :)=[];
        test_y(local_rand)=[];        
    end
    local_flag=local_flag*length(find(training_y==1));
    local_flag=local_flag*length(find(training_y==-1));
    if local_flag==0
        local_flag=1;
    else
        local_flag=0;
    end
end