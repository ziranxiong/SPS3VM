function [x, y, test_x, test_y] = A_read_data(norm)
    global unlabel_size size_training size_data
    load(norm);
    temp = size_training;
    size_training = size_data;
    [x, y, ~, ~] = A_random_select(x, y);
    size_training = temp;
    clear temp
    x = zscore(x);
    [x, y, test_x, test_y] = A_random_select(x, y);
    randomLabelSizeRate = (round(rand(1,1)*9)+1)/100;
    unlabel_size = size_training - randomLabelSizeRate * size_training;
end