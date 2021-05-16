function [time,error,original_x,original_y,extended_x,extended_y,initial_solution,extended_real_y]=select_gamma(data_flag)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [x,y]=read_data(data_flag);
    real_y = y;
    [y]=ModifyData4Semi(y);   
    [original_x,original_y,extended_x,extended_y,extended_real_y]=random_select(x,y,real_y);
    
    time = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Start initial\n');
    initial_solution=initial(original_x,original_y,extended_x,extended_y);
    fprintf('--------------------------------\n');
    %get accuracy
    pred = SemiSVM.predict(initial_solution.x,initial_solution);
    nonLabel_index=initial_solution.label_size+1:length(initial_solution.alpha);
    y1 = extended_real_y(nonLabel_index);
    y2 = pred(nonLabel_index);
    per1 = getAcc(y1,y2);
    error = per1;
end