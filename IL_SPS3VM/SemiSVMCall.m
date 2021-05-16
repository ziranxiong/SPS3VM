function out=SemiSVMCall(original_x,original_y,x,y)
    global C 
    tic;
    label_size=sum(original_y~=0);
    unlabel_size = length(original_y)-label_size;
    os=SMO.InitialSolution(original_x,original_y,C);           
    [PreLowBound,PreUppBound]=SMO.ComputerCCstar(x,y,label_size,os);    
    terminal_flag=0;
    NowLowBound=PreLowBound;
    NowUppBound=PreUppBound;
    tempCount = 0;
    while terminal_flag==0
        objcs=SemiSVM(x,y,NowLowBound,NowUppBound,0);
        PreLowBound=NowLowBound;
        PreUppBound=NowUppBound;
        [NowLowBound,NowUppBound]=SemiSVM.ComputerCCstar(x,y,label_size,objcs);
        terminal_flag=SemiSVM.CheckTerminial(PreLowBound,PreUppBound,NowLowBound,NowUppBound);
        tempCount = tempCount + 1;
    end
    out=objcs;    
    out.local_minimal=SemiSVM.TestMinimal(out,label_size);
    out.label_size = label_size;
    out.unlabel_size = unlabel_size;
    out.time1=toc;
end