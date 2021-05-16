function [time,error,original_x,original_y,extended_x,extended_y,initial_solution,extended_real_y]=main(data_flag)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     [x, y, ~, ~] = A_read_data(norm);
%     real_y = y;
%     [y] = ModifyData4Semi(y);
%     [original_x,original_y,extended_x,extended_y,extended_real_y]=random_select(x,y,real_y);
    
    [x,y]=read_data(data_flag);
    real_y = y;
    [y]=ModifyData4Semi(y);   
    [original_x,original_y,extended_x,extended_y,extended_real_y]=random_select(x,y,real_y);

    
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
    
   % tic;
    [time,error]=IncSemiSVM_algorithm2.Run(original_y,extended_x,extended_y,initial_solution,extended_real_y);

    %get accuracy
%     pred = SemiSVM.predict(out.x,out);
%     nonLabel_index=initial_solution.label_size+1:length(initial_solution.alpha);
%     y1 = extended_real_y(nonLabel_index);
%     y2 = pred(nonLabel_index);
%     per2 = getAcc(y1,y2);
%     
%     out.per1 = per1;
%     out.per2 = per2;
end
