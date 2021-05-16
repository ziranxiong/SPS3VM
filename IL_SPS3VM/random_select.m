function [original_x,original_y,extended_x,extended_y,extended_real_y]=random_select(x,y,real_y)
    global  size_training 
    unlabeled_index=find(y==0);
    labeled_index=[1:length(y)];
    labeled_index(unlabeled_index)=[];
    index_added=[];
    res = [];
    local_size=size_training-length(labeled_index);
    while length(res) ~= local_size
        tmp_res = unique(randi([1,length(unlabeled_index)],local_size,1));
        tmp_res = sort(tmp_res);
        res = [res;tmp_res];
        res = unique(res);
        if length(res) >= local_size
            res = res(randperm(length(res)));
            res = res(1:local_size);
            res = sort(res);
        end
    end
    local_unlabeled_index=unlabeled_index(res);
    final_index=[labeled_index';local_unlabeled_index];
    temp_original_x=x(final_index,:);
    temp_original_y=y(final_index);
    temp_original_real_y = real_y(final_index);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp_nonLabel_index=1:length(temp_original_y);
    temp_Label_index=find(~temp_original_y==0);
    temp_nonLabel_index(temp_Label_index)=[];
    
    %get x0 and y0
    unlabel_len = length(temp_nonLabel_index);
    label_len = length(temp_Label_index);
    x0 = (1/unlabel_len)*sum(temp_original_x(temp_nonLabel_index,:));
    y0 = (1/label_len)*sum(temp_original_y(temp_Label_index));
    
    %add x0,y0 to original_x,original_y
    original_x = [x0;temp_original_x];
    original_y = [y0;temp_original_y];
    original_real_y = [y0;temp_original_real_y];
    
    %get label_index and nonlabel_index
    nonLabel_index=1:length(original_y);
    Label_index=find(~original_y==0);
    nonLabel_index(Label_index)=[];
    
    extended_x=[original_x(Label_index,:);original_x(nonLabel_index,:);original_x(nonLabel_index,:)];
    extended_y=[original_y(Label_index);-ones(length(nonLabel_index),1);ones(length(nonLabel_index),1)];
    extended_real_y = [original_real_y(Label_index);original_real_y(nonLabel_index);original_real_y(nonLabel_index)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end