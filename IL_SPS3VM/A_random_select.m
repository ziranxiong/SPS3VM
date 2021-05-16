function [training_x, training_y, test_x, test_y]=A_random_select(x, y)
% 随机选取数据
% 将数据集划分为 测试集(test) & 训练集(training)
% 测试集：test_x, test_y
% 训练集：training_x, training_y
% 使用 size_training 控制训练集的大小
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
        % 随机产生一个 1 * size_training-len1 的矩阵，矩阵中的数据在[1, length-len1]之间
        local_rand= randi([1,local_length-len1],1,size_training-len1);
        % 去掉矩阵中的重复数据
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