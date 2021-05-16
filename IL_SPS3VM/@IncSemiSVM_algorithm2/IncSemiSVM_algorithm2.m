classdef IncSemiSVM_algorithm2  
    properties
        KKT=0;
    end
    methods(Static = true)
        function obj=IncSemiSVM_algorithm2(a)
            obj.KKT=a;
        end
        function [time,error]=Run(original_y,extended_x,extended_y,initial_solution,extended_real_y)
            global fake_zero C_end C index_Q flag C_star s
            index_ox = [1:length(original_y)];
%             obj.o_x=initial_solution.x(index_ox,:);
%             obj.o_y=original_y;
%             obj.ex_x=extended_x;
           % obj.ex_y=extended_y;
            flag=0;
            svm_struct=initial_solution;
            
%             preSumX = svm_struct.x(1,:)*svm_struct.unlabel_size;
%             preSumY = svm_struct.y(1)*svm_struct.label_size;
                
            % remoeve x_0 y_0
            fprintf('step 1: remoeve x_0 y_0\n');
            remove_y = svm_struct.y(1);
            g=IncSemiSVM_algorithm2.GetG(svm_struct.x,svm_struct.y,svm_struct);    %get gc   
            [svm_struct,index_A,direction_A]=IncSemiSVM_algorithm2.InitialA2(svm_struct,g(1),1,remove_y); 
            [Qc]=IncSemiSVM_algorithm2.GetQc(svm_struct,svm_struct.x,index_A,direction_A);
            while svm_struct.alpha(1) > fake_zero || svm_struct.alpha(1) < -fake_zero
                [beta]=IncSemiSVM_algorithm2.GetBeta(extended_x,extended_y,svm_struct,Qc);
                [gama] = IncSemiSVM_algorithm2.GetGama(extended_x,extended_y,svm_struct,beta,direction_A,index_A);
                [max_delta,index,index_class]=IncSemiSVM_algorithm2.GetMaxChange2...
                    (svm_struct,g,beta,gama,extended_y,direction_A,index_A,original_y);
%                 if(old_max_delta~=-9999)
%                    if(old_max_delta~=max_delta)
%                        change_zeta=change_zeta+1;
%                    end
%                 end
%                 old_max_delta=max_delta;
                [flag,g,svm_struct,Qc,direction_A,index_A]=IncSemiSVM_algorithm2.Update...
                    (svm_struct,g,beta,gama,max_delta,index,extended_x,extended_y,direction_A,index_A,index_class);               
%                 if flag~=0
%                     break
%                 end
%                 steps=steps+1;
            end
            x_0 = svm_struct.x(1,:);
            y_0 = svm_struct.y(1);
            svm_struct.alpha(1) = [];
            prePms0 = svm_struct.pms(1,:);
            svm_struct.pms(1,:) = [];
            svm_struct.x(1,:) = [];
            svm_struct.y(1) = [];
            g(1)=[];
            svm_struct.label_size = svm_struct.label_size - 1;
            svm_struct.gi=g;
            
            initLowIn = svm_struct.pms(:,1);
            initUbIn = svm_struct.pms(:,2);
            svm_struct.index_M=find((svm_struct.alpha>initLowIn+fake_zero) & (svm_struct.alpha<initUbIn-fake_zero)); 
            svm_struct.index_E=find( (svm_struct.alpha>=initUbIn-fake_zero));
            svm_struct.index_O=find( (svm_struct.alpha<=initLowIn+fake_zero));
            % SPS3VM
%             fprintf('step 2: SPS3VM\n');
%                  error=[];
%                  time=[];
%             index_Q=[svm_struct.index_E(svm_struct.alpha(svm_struct.index_E)~=0 & svm_struct.alpha(svm_struct.index_E)~=C);svm_struct.index_O(svm_struct.alpha(svm_struct.index_O)~=0)];
%             d=C_end-C_star;
%             [SPQC]=IncSemiSVM_algorithm2.SPGetQc(svm_struct,d);
%             while (C_star<C_end)
%                 tic;
%                 if isempty(index_Q)
%                     fprintf('Q is empty\n');
%                     break;
%                 end
%                 [beta]=IncSemiSVM_algorithm2.SPGetBeta(svm_struct.x,svm_struct.y,svm_struct,SPQC);
%                 [gama]=IncSemiSVM_algorithm2.SPGetGama(svm_struct.x,svm_struct.y,svm_struct,beta,d);
%                 [max_delta,index,index_class]=IncSemiSVM_algorithm2.SPGetMaxChange...
%                     (svm_struct,g,beta,gama,svm_struct.y,d,svm_struct.y);
%                 [g,svm_struct,SPQC,d]=IncSemiSVM_algorithm2.SPUpdate...
%                     (svm_struct,g,beta,gama,max_delta,index,svm_struct.x,svm_struct.y,d,index_class);
%                 t1=toc;
%                 test_svm=svm_struct; 
%                 if flag==1
%                     break;
%                 end

            Q=SMO.Kernel(svm_struct.x,svm_struct.x).*(svm_struct.y*(svm_struct.y)');
            local_Q=[svm_struct.y Q];
            local_alpha=[svm_struct.b;svm_struct.alpha];
            local_f=local_Q*local_alpha;
            index_Q=find(local_f>s & local_f<0);
        
            while length(index_Q)~=0
                MM=local_f(index_Q);
                [p,q]=sort(MM);
                index=index_Q(q(1));
                s=local_f(index)+1e-10;
                svm_struct.alpha(index) = [];
                prePms0 = svm_struct.pms(index,:);
                svm_struct.pms(index,:) = [];
                svm_struct.x(index,:) = [];
                svm_struct.y(index) = [];
                g(index)=[];
                %svm_struct.label_size = svm_struct.label_size - 1;
                svm_struct.gi=g;

                initLowIn = svm_struct.pms(:,1);
                initUbIn = svm_struct.pms(:,2);
                svm_struct.index_M=find((svm_struct.alpha>initLowIn+fake_zero) & (svm_struct.alpha<initUbIn-fake_zero)); 
                svm_struct.index_E=find( (svm_struct.alpha>=initUbIn-fake_zero));
                svm_struct.index_O=find( (svm_struct.alpha<=initLowIn+fake_zero));
%                 % add x_0 y_0
%                 %fprintf('step 3: add x_0 y_0 \n')
%                 test_svm.label_size = test_svm.label_size + 1;
%                 %init a0,x0,y0
   

out = IncSemiSVM.Run(original_y,extended_x,extended_y,added_x,added_y,initial_solution);

alpha0 = 0;
%                 %add sample 0

                test_svm=svm_struct;
                test_svm.x = [test_svm.x;x_0];
                test_svm.y = [test_svm.y;y_0];
                test_svm.alpha = [test_svm.alpha;0];
               
%             Q=SMO.Kernel(svm_struct.x,svm_struct.x).*(svm_struct.y*(svm_struct.y)');
%             local_Q=[svm_struct.y Q];
%             local_alpha=[svm_struct.b;svm_struct.alpha];
%             local_f=local_Q*local_alpha;
%   
%             mu=C_star*(local_f<s);
%             low_tmp=(-mu).*(y>0)+(mu-C).*(y<0);
%             upp_tmp=(C-mu).*(y>0)+(mu).*(y<0);
                


%             LowBoud2=zeros(test_svm.label_size,1);
%             UppBound2=C*ones(test_svm.label_size,1);
            Q=SMO.Kernel(test_svm.x,test_svm.x).*(test_svm.y*(test_svm.y)');
            local_Q=[test_svm.y Q];
            local_alpha=[test_svm.b;test_svm.alpha];
            local_f=local_Q*local_alpha;
            mu=C_star*(local_f<s);
            low_tmp=(-mu).*(test_svm.y>0)+(mu-C).*(test_svm.y<0);
            upp_tmp=(C-mu).*(test_svm.y>0)+(mu).*(test_svm.y<0);
%             LowBoud2=[LowBoud2;low_tmp];
%             UppBound2=[UppBound2;upp_tmp];


%                 prePms0(1) = low_tmp;
%                 prePms0(2) = upp_tmp;
                test_svm.pms=[low_tmp upp_tmp];
                
                initLowIn = test_svm.pms(:,1);
                initUbIn = test_svm.pms(:,2);
                test_svm.index_M=find((test_svm.alpha>initLowIn+fake_zero) & (test_svm.alpha<initUbIn-fake_zero));
                test_svm.index_E=find( (test_svm.alpha>=initUbIn-fake_zero));
                test_svm.index_O=find( (test_svm.alpha<=initLowIn+fake_zero));
                
                
                g1=IncSemiSVM_algorithm2.GetG(test_svm.x,test_svm.y,test_svm);    %get gc
               
                
                [test_svm,index_A,direction_A]=IncSemiSVM_algorithm2.InitialA3(test_svm,g(end),length(g),y_0);
               
                
                [Qc]=IncSemiSVM_algorithm2.GetQc(test_svm,test_svm.x,index_A,direction_A);
%                 if index_A>1
%                     fprintf('add x0 y0 error\n');
%                 end
tic;
                while length(index_A)>0
                    [beta_new]=IncSemiSVM_algorithm2.GetBeta(extended_x,extended_y,test_svm,Qc);
                    [gama_new] = IncSemiSVM_algorithm2.GetGama(extended_x,extended_y,test_svm,beta_new,direction_A,index_A);
                    [max_delta,index,index_class]=IncSemiSVM_algorithm2.GetMaxChange...
                        (test_svm,g1,beta_new,gama_new,test_svm.y,direction_A,index_A,original_y);
                    [flag,g1,test_svm,Qc,direction_A,index_A]=IncSemiSVM_algorithm2.Update...
                        (test_svm,g1,beta_new,gama_new,max_delta,index,extended_x,extended_y,direction_A,index_A,index_class);
%                     if flag~=0
%                         break
%                     end 
                end
                t1=toc;
                Q=SMO.Kernel(svm_struct.x,svm_struct.x).*(svm_struct.y*(svm_struct.y)');
                local_Q=[svm_struct.y Q];
                local_alpha=[svm_struct.b;svm_struct.alpha];
                local_f=local_Q*local_alpha;
                index_Q=find(local_f>s) && find(local_f<0);
            end
                pred = SemiSVM.predict(test_svm.x,test_svm);
                nonLabel_index=test_svm.label_size+1:length(test_svm.alpha);
                y1 = extended_real_y(nonLabel_index);
                y2 = pred(nonLabel_index);
                per1 = getAcc(y1,y2);
                error=[error;per1];
                time=[time;t1];
        end

            %using IL-S3VM again end -------------increment alpha0
%             [svm_struct,label_size]=IncSemiSVM_algorithm2.ModifyDataAfterInc(original_y,added_y,svm_struct);
%             out_KKT1=SemiSVM.TestKKT(svm_struct,svm_struct.pms(:,1),svm_struct.pms(:,2));
%             svm_struct.KKT=out_KKT1;
%             svm_struct.local_minimal=SemiSVM.TestMinimal(svm_struct,label_size);
%             svm_struct.steps=steps;
%             svm_struct.change_zeta=change_zeta;
%             svm_struct.flag=flag;
%             out=svm_struct;
        
        function   [flag,out_g,out_svm_struct,out_Qc,direction_A,index_A]=Update...
                (svm_struct,g,beta,gama,max_delta,index,x,y,direction_A,index_A,index_class)
            [flag,out_svm_struct,direction_A,index_A]=IncSemiSVM_algorithm2.UpdataSvmstruct(svm_struct,beta,max_delta,index,y,direction_A,index_A,index_class);
            [out_g]=IncSemiSVM_algorithm2.UpdataG(g,svm_struct,gama,max_delta);
            [out_Qc]=IncSemiSVM_algorithm2.GetQc(out_svm_struct,x,index_A,direction_A);
        end
        function [out_g]=UpdataG(g,svm_struct,gama,max_delta)
            other_set=[1:length(g)];
            other_set(svm_struct.index_M)=[];
            g(other_set)=g(other_set)+gama(other_set)*max_delta;
            out_g=g;
        end
        function  [flag,out_svm_struct,direction_A,index_A] = UpdataSvmstruct...
                (svm_struct,beta,max_delta,index,y,direction_A,index_A,index_class)
            global  fake_zero size_training C_star
            flag=0;
            local_beta=beta;
            local_beta_length=length(local_beta);
            beta_b=beta(1);
            beta_s=beta(2:local_beta_length);
            alpha=svm_struct.alpha;
            support_set=svm_struct.index_M;
            b=svm_struct.b;
            alpha(support_set)=alpha(support_set)+beta_s*max_delta;
            alpha(index_A) = alpha(index_A)+direction_A*max_delta;
            b=b+beta_b*max_delta;
            svm_struct.b=b;
            svm_struct.alpha=alpha;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% details for index_class
            %  1 : alpha moves from M set to O set
            %  2 :  alpha moves from M set to E set
            %  3 :  g moves from O set to M set
            %  4 :  g moves from E set to M set
            %  5 :  a sample moves from A to M
            %  6 :  alpha moves from A set to O set
            %  7 :  alpha moves from A set to E set
            %  8 :  sample moves from O set to A set
            %  9 :  g moves from E set to A set
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%set
            switch  index_class
                case 1
                    tmp_index = find(svm_struct.index_M == index);
                    svm_struct.index_M(tmp_index)=[];
                    svm_struct.index_O = sort([svm_struct.index_O;index]);
                case 2
                    tmp_index = find(svm_struct.index_M == index);
                    svm_struct.index_M(tmp_index)=[];
                    svm_struct.index_E = sort([svm_struct.index_E;index]);
                case 3
                    tmp_index = find(svm_struct.index_O == index);
                    svm_struct.index_O(tmp_index)=[];
                    svm_struct.index_M = sort([svm_struct.index_M;index]);
                case 4
                    tmp_index = find(svm_struct.index_E == index);
                    svm_struct.index_E(tmp_index)=[];
                    svm_struct.index_M = sort([svm_struct.index_M;index]);
                case 5
                    tmp_index = find(index_A == index);
                    index_A(tmp_index)=[];
                    direction_A(tmp_index,:) = [];
                    svm_struct.index_M = sort([svm_struct.index_M;index]);
                case 6
                    tmp_index = find(index_A == index);
                    index_A(tmp_index)=[];
                    direction_A(tmp_index,:) = [];
                    svm_struct.index_O = sort([svm_struct.index_O;index]);
                case 7
                    tmp_index = find(index_A == index);
                    index_A(tmp_index)=[];
                    direction_A(tmp_index,:) = [];
                    svm_struct.index_E = sort([svm_struct.index_E;index]);
                case 8
                    tmp_index = find(svm_struct.index_O == index);
                    svm_struct.index_O(tmp_index)=[];
                    index_A = [index_A;index];                   
                    tmp1 = svm_struct.pms(index,2)*(svm_struct.y(index)>0)+svm_struct.pms(index,1)*(svm_struct.y(index)<0);
                    local_direction = tmp1-svm_struct.alpha(index);
                    direction_A = [direction_A;local_direction];
                case 9
                    tmp_index = find(svm_struct.index_E == index);
                    svm_struct.index_E(tmp_index)=[];
                    index_A = [index_A; index];
                    
                    tmp1 = svm_struct.pms(index,2)*(svm_struct.y(index)>0)+svm_struct.pms(index,1)*(svm_struct.y(index)<0);
                    
                    local_direction = tmp1-svm_struct.alpha(index);
                    direction_A = [direction_A;local_direction];
                case 10  
                    tmp1 = svm_struct.pms(index,2)*(svm_struct.y(index)>0)+svm_struct.pms(index,1)*(svm_struct.y(index)<0);
                    local_direction = tmp1-svm_struct.alpha(index);
                    tmp2=find(index_A==index);
                    direction_A(tmp2) = local_direction;
                otherwise
                    %error('has a error for class_toA');
                    flag=1;
            end
            out_svm_struct=svm_struct;
        end
 
        function  [max_delta_zeta,index,index_class]=GetMaxChange...
                (svm_struct,g,beta,gama,y,direction_A,index_A,original_y)
            [max_change_alpha,index_max_change_alpha,class_alpha]=IncSemiSVM_algorithm2.FindChangeAlpha(svm_struct,beta);
            [max_change_g,index_max_change_g,class_g]=IncSemiSVM_algorithm2.FindMaxChangeG(svm_struct,g,gama);
            [max_change_A,index_max_change_A,class_A]=IncSemiSVM_algorithm2.FindChangeA...
                (svm_struct,y,beta, direction_A,index_A,g,gama);
            [max_change_toA,index_max_change_toA,class_toA]=IncSemiSVM_algorithm2.FindChangeToA...
                (svm_struct,y, index_A,g,gama,original_y);
            [max_delta_zeta,index,index_class]=IncSemiSVM_algorithm2.FindMaxChange...
                (max_change_alpha, index_max_change_alpha,max_change_g,...
                index_max_change_g,max_change_A, index_max_change_A, max_change_toA, ...
                index_max_change_toA,class_alpha,class_g,class_A,class_toA);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% details for class_toA
            %  1 : alpha moves from M set to O set
            %  2 :  alpha moves from M set to E set
            %  3 :  g moves from O set to M set
            %  4 :  g moves from E set to M set
            %  5 :  a sample moves from A to M
            %  6 :  alpha moves from A set to O set
            %  7 :  alpha moves from A set to E set
            %  8 :  sample moves from O set to A set
            %  9 :  g moves from E set to A set
            %  10:  g moves from A set to A set
        end

        function  [max_delta_zeta,index,index_class]=GetMaxChange2...
                (svm_struct,g,beta,gama,y,direction_A,index_A,original_y)
            [max_change_alpha,index_max_change_alpha,class_alpha]=IncSemiSVM_algorithm2.FindChangeAlpha(svm_struct,beta);
            [max_change_g,index_max_change_g,class_g]=IncSemiSVM_algorithm2.FindMaxChangeG(svm_struct,g,gama);
            [max_change_A,index_max_change_A,class_A]=IncSemiSVM_algorithm2.FindChangeA2...
                (svm_struct,y,beta, direction_A,index_A,g,gama);
            [max_change_toA,index_max_change_toA,class_toA]=IncSemiSVM_algorithm2.FindChangeToA...
                (svm_struct,y, index_A,g,gama,original_y);
            [max_delta_zeta,index,index_class]=IncSemiSVM_algorithm2.FindMaxChange...
                (max_change_alpha, index_max_change_alpha,max_change_g,...
                index_max_change_g,max_change_A, index_max_change_A, max_change_toA, ...
                index_max_change_toA,class_alpha,class_g,class_A,class_toA);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% details for class_toA
            %  1 : alpha moves from M set to O set
            %  2 :  alpha moves from M set to E set
            %  3 :  g moves from O set to M set
            %  4 :  g moves from E set to M set
            %  5 :  a sample moves from A to M
            %  6 :  alpha moves from A set to O set
            %  7 :  alpha moves from A set to E set
            %  8 :  sample moves from O set to A set
            %  9 :  g moves from E set to A set
            %  10:  g moves from A set to A set
        end
        function    [updata_max_change, updata_index,index_class]=FindMaxChange...
                (max_change_alpha, index_max_change_alpha,max_change_g,...
                index_max_change_g,max_change_A, index_max_change_A, max_change_toA, index_max_change_toA, ...
                class_alpha,class_g,class_A,class_toA)
            updata_max_change=inf;
            updata_index=0;
            index_class=inf;
            if max_change_alpha <= updata_max_change
                updata_max_change=max_change_alpha;
                updata_index=index_max_change_alpha;
                index_class=class_alpha;
            end
            if max_change_g  <= updata_max_change 
                updata_max_change=max_change_g;
                updata_index=index_max_change_g;
                index_class=class_g;
            end
            if max_change_A<=updata_max_change
                updata_max_change=max_change_A;
                updata_index=index_max_change_A;
                index_class=class_A;
            end
            if max_change_toA <=updata_max_change
                updata_max_change=max_change_toA;
                updata_index=index_max_change_toA;
                index_class=class_toA;
            end
        end
        function [max_change_A, index_max_change_A, class_A] = FindChangeA...
                (svm_struct,y,beta,direction_A,index_A,g,gama)
            % condition 4 in the paper
            global fake_zero
            class_A=inf;
            local_g=g(index_A);
            local_gama=gama(index_A);
            local_index_A = index_A;
            local_direction_A =direction_A;
            max_change_A=inf;
            index_max_change_A=0;
            local_index=find(local_gama<fake_zero & local_gama>-fake_zero);
            local_gama(local_index)=[];
            local_index_A(local_index)=[];
            local_direction_A(local_index)=[];
            local_g(local_index)=[];
            local_index3=find(local_g.*local_gama>=0);
            local_g(local_index3)=[];
            local_gama(local_index3)=[];
            local_index_A(local_index3)=[];
            local_direction_A(local_index3)=[];
            local_change_g=-local_g./local_gama;
            low_bound=svm_struct.pms(local_index_A,1);
            upp_bound=svm_struct.pms(local_index_A,2);
            if length(local_index_A)>0
                tmp5 = svm_struct.alpha(local_index_A)+local_direction_A.*local_change_g;
                local_index4=find((tmp5<low_bound) | (tmp5>upp_bound));
                local_index_A(local_index4)=[];
                local_change_g(local_index4)=[];
            end
            [local_change local_index]=min(local_change_g);
            if length(local_change_g)==0
                max_change_A=inf;
            else
                max_change_A=local_change;
                index_max_change_A=local_index_A(local_index);
                class_A = 5; %%% a sample moves from A to M
            end
            low_bound=svm_struct.pms(index_A,1);
            upp_bound=svm_struct.pms(index_A,2);
            local_index_A=index_A;
            local_direction_A=direction_A;
            low_terminal = low_bound;
            if length(local_index_A)>0
                local_change_al = (low_terminal-svm_struct.alpha(local_index_A))./local_direction_A;
                local_index3=find(local_change_al<0);
                local_change_al(local_index3)=[];
                local_index_A(local_index3)=[];
            end
            if length(local_index_A)>0
                tmpg6 = g(local_index_A)+gama(local_index_A).*local_change_al;
                local_index4=find(tmpg6 <0);
                local_index_A(local_index4)=[];
                local_change_al(local_index4)=[];
            end
            [min_alpha_change,  min_alpha_change_index]=min(local_change_al);
            if min_alpha_change < max_change_A
                max_change_A = min_alpha_change;
                index_max_change_A=local_index_A(min_alpha_change_index);
                class_A=6;  % alpha moves from A set to O set
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5               
            local_direction_A=direction_A;
            upp_terminal = upp_bound;
            local_index_A=index_A;
            if length(local_index_A)>0
                local_change_al = (upp_terminal-svm_struct.alpha(local_index_A))./local_direction_A;
                local_index3=find(local_change_al<0);
                local_change_al(local_index3)=[];
                local_index_A(local_index3)=[];
            end
            if length(local_index_A)>0
                tmpg6 = g(local_index_A)+gama(local_index_A).*local_change_al;
                local_index4=find(tmpg6 > 0);
                local_index_A(local_index4)=[];
                local_change_al(local_index4)=[];
            end
            [min_alpha_change  min_alpha_change_index]=min(local_change_al);
            if min_alpha_change < max_change_A
                max_change_A = min_alpha_change;
                index_max_change_A=local_index_A(min_alpha_change_index);
                class_A=7;  % alpha moves from A set to E set
            end
        end
        function [max_change_A, index_max_change_A, class_A] = FindChangeA2...
                (svm_struct,y,beta,direction_A,index_A,g,gama)
            % condition 4 in the paper
            global fake_zero
            class_A=inf;
            local_g=g(index_A);
            local_gama=gama(index_A);
            local_index_A = index_A;
            local_direction_A =direction_A;
            max_change_A=inf;
            index_max_change_A=0;
            local_index=find(local_gama<fake_zero & local_gama>-fake_zero);
            local_gama(local_index)=[];
            local_index_A(local_index)=[];
            local_direction_A(local_index)=[];
            local_g(local_index)=[];
            low_bound=svm_struct.pms(local_index_A,1);
            upp_bound=svm_struct.pms(local_index_A,2);
            max_change_A = inf;
            if length(local_index_A)>0
                local_change_g = (low_bound - svm_struct.alpha(local_index_A))./local_direction_A
                max_change_A = local_change_g
            end
            index_max_change_A = 1
        end
        function [max_change_toA,index_max_change_toA, class_toA]=FindChangeToA...
                (svm_struct,y, index_A,g,gama,original_y)
            % condition 5 in the paper
            global fake_zero
            class_toA = inf;
            len_labeled = length(find(original_y~=0));
            index_labeled = [1:len_labeled]';
            tmp1=unique([index_labeled;svm_struct.index_M;]);
            other_set=[1:length(y)];
            other_set(tmp1)=[];
            gama=gama(other_set);
            local_g=g(other_set);
            max_change_toA=inf;
            index_max_change_toA=0;
            local_index=find(gama<fake_zero & gama>-fake_zero);
            gama(local_index)=[];
            other_set(local_index)=[];
            local_g(local_index)=[];
            local_g=local_g.*y(other_set)+1;
            gama=gama.*y(other_set);            
            local_index3=find(local_g.*gama>=0);
            local_g(local_index3)=[];
            gama(local_index3)=[];
            other_set(local_index3)=[];
            local_change_g = -local_g./gama;
            [local_change, local_index]=min(local_change_g);
            if isempty(local_change_g)
                max_change_toA=inf;
            else
                max_change_toA=local_change;
                index_max_change_toA=other_set(local_index);
                if length(find(index_A==index_max_change_toA))>0
                    class_toA = 10;    % g moves from A set to A set
                else
                    if g(index_max_change_toA) < 0
                        class_toA = 9;  % g moves from E set to A set
                    else
                        class_toA = 8;   % sample moves from O set to A set
                    end
                end
            end
        end
        function [max_change_alpha, index_max_change_alpha,class_alpha] = FindChangeAlpha(svm_struct,beta)
            global fake_zero
            index_max_change_alpha=0;
            class_alpha=inf;
            support_set=svm_struct.index_M;
            alpha=svm_struct.alpha;
            local_alpha=alpha(support_set);
            local_beta_length=length(beta);
            local_beta=beta(2:local_beta_length);
            local_index=find(local_beta<fake_zero & local_beta>-fake_zero);
            local_alpha(local_index)=[];
            support_set(local_index)=[];
            local_beta(local_index)=[];
            local_C=svm_struct.pms(support_set,1).*(local_beta<0)+svm_struct.pms(support_set,2).*(local_beta>0);
            local_change_al=(local_C-local_alpha)./local_beta;
            [min_alpha_change, min_alpha_change_index]=min(local_change_al);
            if length(local_change_al)==0
                max_change_alpha=inf;
            else
                max_change_alpha=min_alpha_change;
                index_max_change_alpha=support_set(min_alpha_change_index);
                if local_beta(min_alpha_change_index)>0
                    class_alpha=2;  % alpha moves from M set to E set
                else
                    class_alpha=1;  % alpha moves from M set to O set
                end
            end
        end
        function [max_change_g,index_max_change_g,class_g] = FindMaxChangeG(svm_struct,g,gama)
            global fake_zero
            class_g = inf;
            other_set=[1:length(g)];
            other_set(svm_struct.index_M)=[];
%             other_set(end)=[];  %% delete the candidate sample from the other set
            gama=gama(other_set);        %% delete the candidate sample from the gama
            local_g=g(other_set);
            max_change_g=inf;
            index_max_change_g=0;
            local_index=find(gama<fake_zero & gama>-fake_zero);
            gama(local_index)=[];
            other_set(local_index)=[];
            local_g(local_index)=[];
            local_index2=find(local_g<fake_zero & local_g>-fake_zero);
            local_g(local_index2)=[];
            gama(local_index2)=[];
            other_set(local_index2)=[];
            local_index3=find(local_g.*gama>=0);
            local_g(local_index3)=[];
            gama(local_index3)=[];
            other_set(local_index3)=[];
            local_change_g=-local_g./gama;
            [local_change local_index]=min(local_change_g);
            if length(local_change_g)==0
                max_change_g=inf;
            else
                max_change_g=local_change;
                index_max_change_g=other_set(local_index);
                if gama(local_index)>0
                    class_g=4;  % g moves from E set to M set
                else
                    class_g=3;  % g moves from O set to M set
                end
            end
        end
        function  g=GetG(x,y,svm_struct)
            Q=SemiSVM.Kernel(x,svm_struct.x);
            wx=Q*svm_struct.alpha;
            g=wx + svm_struct.b-y;
        end
        function [svm_struct,index_A,direction_A] = InitialA(svm_struct,gc,index_added,added_y)
            index_A=[];
            direction_A=[];
            alpha_c=svm_struct.alpha(end);
            pms=svm_struct.pms(end,:);
            global fake_zero
            if gc<-fake_zero & alpha_c>pms(2)-fake_zero
                svm_struct.index_E=[svm_struct.index_E;index_added];
            else
                if gc>fake_zero & alpha_c<pms(1)+fake_zero
                    svm_struct.index_O=[svm_struct.index_O;index_added];
                else
                    index_A=[index_added];
                    index_M=svm_struct.index_M;
                    length_A=length(index_A);
                    y_A=svm_struct.y(index_A);
                    alpha=svm_struct.alpha;
                    pms=svm_struct.pms;
                    if added_y==0
                        tmp1 = (y_A>0).*pms(index_A,2)+(y_A<0).*pms(index_A,1);
                    else
                        tmp1 = pms(index_A,2);
                    end
                    direction_A=tmp1-alpha(index_A);
                end
            end            
        end
        function [svm_struct,index_A,direction_A] = InitialA2(svm_struct,gc,index_remove,remove_y)
            index_A=[];
            direction_A=[];
            alpha_c=svm_struct.alpha(1);
            pms=svm_struct.pms(1,:);
            global fake_zero
            %if gc<-fake_zero & alpha_c>pms(2)-fake_zero
            %    svm_struct.index_E=[svm_struct.index_E;index_remove];
            %else
            %    if gc>fake_zero & alpha_c<pms(1)+fake_zero
            %        svm_struct.index_O=[svm_struct.index_O;index_remove];
            %    else
                    index_A=[index_remove];
                    index_M=svm_struct.index_M;
                    length_A=length(index_A);
                    y_A=svm_struct.y(index_A);
                    alpha=svm_struct.alpha;
                    pms=svm_struct.pms;
                    if remove_y==0
                        tmp1 = (y_A>0).*pms(index_A,1)+(y_A<0).*pms(index_A,2);
                    else
                        tmp1 = pms(index_A,1);
                    end
                    direction_A=tmp1-alpha(index_A);
             %   end
            %end            
        end
        function [svm_struct,index_A,direction_A] = InitialA3(svm_struct,gc,index_added,added_y)
            index_A=[];
            direction_A=[];
            alpha_c=svm_struct.alpha(1);
            pms=svm_struct.pms(1,:);
            global fake_zero
            if gc<-fake_zero & alpha_c>pms(2)-fake_zero
                svm_struct.index_E=[svm_struct.index_E;index_added];
            else
                if gc>fake_zero & alpha_c<pms(1)+fake_zero
                    svm_struct.index_O=[svm_struct.index_O;index_added];
                else
                    index_A=[index_added];
                    index_M=svm_struct.index_M;
                    length_A=length(index_A);
                    y_A=svm_struct.y(index_A);
                    alpha=svm_struct.alpha;
                    pms=svm_struct.pms;
                    if added_y==0
                        tmp1 = (y_A>0).*pms(index_A,2)+(y_A<0).*pms(index_A,1);
                    else
                        tmp1 = pms(index_A,2);
                    end
                    direction_A=tmp1-alpha(index_A);
                end
            end            
        end
        function [Qc]=GetQc(svm_struct,x,index_A,direction_A)
            Qc=[];
            local_direction=direction_A;
            if length(index_A)>0
                index_M=svm_struct.index_M;
                length_A=length(index_A);
%                 y_A=svm_struct.y(index_A);
%                 alpha=svm_struct.alpha;
%                 pms=svm_struct.pms;
                tmp1=ones(1,length_A);
                x_M=x(index_M,:);
                x_A=x(index_A,:);
                H_MA=SemiSVM.Kernel(x_M,x_A);
                tmp_H=[tmp1;H_MA];
                Qc=-tmp_H*local_direction;
            end
        end
        function    [gama]=GetGama(x,y,svm_struct,beta,direction_A,index_A)
%             other_set=[1:length(y)];
%             other_set(svm_struct.index_M)=[];
%             length_other_set=length(other_set);
            support_set=svm_struct.index_M;
            index_V=svm_struct.index_V;
            local_x=x;
            local_sx=x(support_set,:);
            local_Ax=x(index_A,:);
            local_Vx=x(index_V,:);
            local_Mbeta=SemiSVM.Kernel(local_x,local_sx);
            local_Vbeta=SemiSVM.Kernel(local_x,local_Vx);
            local_qA=SemiSVM.Kernel(local_x,local_Ax);
            local_qbeta=[ones(length(y),1)  local_Mbeta local_Vbeta];
            gama=local_qbeta*beta+local_qA*direction_A;
        end
        function    beta=GetBeta(x,y,svm_struct,Qc)
            Q=IncSemiSVM_algorithm2.FindQR(x,y,svm_struct);
            if det(Q)<1e-42
                m=length(Q);
                n=eye(m);
                Q=Q+1e-7*n;
            end
            beta=linsolve(Q,Qc);
            
            %beta=linsolve(Q,Qc);
        end
        function [svm_struct,label_size]=ModifyDataAfterInc(original_y,added_y,svm_struct)
            label_size=sum(original_y~=0);
            if added_y~=0
                tmp_x=svm_struct.x(end,:);
                svm_struct.x(end,:)=[];
                tmp_y=svm_struct.y(end,:);
                svm_struct.y(end)=[];
                svm_struct.x=[svm_struct.x([1:label_size],:);tmp_x;svm_struct.x([label_size+1:end],:)];
                svm_struct.y=[svm_struct.y([1:label_size]);tmp_y;svm_struct.y([label_size+1:end])];
                svm_struct.alpha=[svm_struct.alpha([1:label_size]);svm_struct.alpha(end);svm_struct.alpha([label_size+1:end-1])];
                svm_struct.pms=[svm_struct.pms([1:label_size],:);svm_struct.pms(end,:);svm_struct.pms([label_size+1:end-1],:)];
                tmp1=find(svm_struct.index_M==length(svm_struct.y));
                if sum(svm_struct.index_M>label_size)>0
                    svm_struct.index_M(svm_struct.index_M>label_size)=svm_struct.index_M(svm_struct.index_M>label_size)+1;
                end
                if length(tmp1)>0
                    svm_struct.index_M(tmp1)=label_size+1;
                end
                tmp1=find(svm_struct.index_O==length(svm_struct.y));
                if sum(svm_struct.index_O>label_size)>0
                    svm_struct.index_O(svm_struct.index_O>label_size)=svm_struct.index_O(svm_struct.index_O>label_size)+1;
                end
                if length(tmp1)>0
                    svm_struct.index_O(tmp1)=label_size+1;
                end
                tmp1=find(svm_struct.index_E==length(svm_struct.y));
                if sum(svm_struct.index_E>label_size)>0
                    svm_struct.index_E(svm_struct.index_E>label_size)=svm_struct.index_E(svm_struct.index_E>label_size)+1;
                end
                if length(tmp1)>0
                    svm_struct.index_E(tmp1)=label_size+1;
                end
                label_size=label_size+1;
            end
        end
        
        function Q=FindQR(x,y,svm_struct)
            global  NumSingular fake_zero
            index_M=svm_struct.index_M;
            length_M=length(index_M);
            x_M=x(index_M,:);
            Q=SemiSVM.Kernel(x_M,x_M);
            tmp1=ones(length_M,1);
            tmp1p=[0;tmp1];
            Q=[tmp1';Q];
            Q=[tmp1p  Q];
            tmp=det(Q);
%             if tmp<fake_zero & tmp>-fake_zero
%                 NumSingular=NumSingular+1;
%             end
        end
        function Q=FindQR2(x,y,svm_struct)
            global  NumSingular fake_zero
            index_M=svm_struct.index_M;
            index_V=svm_struct.index_V;
            length_M=length(index_M);
            length_V=length(index_V);
            x_M=x(index_M,:);
            x_V=x(index_V,:);
            H_MV=SemiSVM.Kernel(x_M,x_V);
            H_MM=SemiSVM.Kernel(x_M,x_M);
            H_1M=ones(length_M,1);
            H_1V=ones(length_V,1);
            tmp1p=[0;H_1M];
            tem2p=[H_1M';H_MM];
            tem3p=[H_1V';H_MV];
            Q=[tmp1p tem2p tem3p];
%            tmp=det(Q);
%             if tmp<fake_zero & tmp>-fake_zero
%                 NumSingular=NumSingular+1;
%             end
        end
        function out_KKT=TestKKT(svm_struct,g)
            global fake_zero2
            zero=sum(svm_struct.alpha);
            out_KKT=0;
            if zero<fake_zero2 && zero>-fake_zero2
                index=sort([svm_struct.index_M;svm_struct.index_O;svm_struct.index_E]);
                alpha=svm_struct.alpha(index);
                local_g=g(index);
                C_low=svm_struct.pms(index,1);
                C_up=svm_struct.pms(index,2);
                out_KKT=SemiSVM.SubKKT(alpha,local_g,C_low,C_up,index);
            end
        end
        
        %SPS3VM  about function      
        function [Qc]=SPGetQc(svm_struct,d)
            global  index_Q
            x=svm_struct.x;
            index_M=svm_struct.index_M;
            x_M=x(index_M,:);
            x_Q=x(index_Q,:);
            H1=ones(1,length(index_Q));
            H2=SemiSVM.Kernel(x_M,x_Q);
            y_Qd=svm_struct.y(index_Q)*d;
            tmp_H=[H1;H2]*y_Qd;
            Qc=-tmp_H;
        end
        function [gama]=SPGetGama(x,y,svm_struct,beta,d)
            global index_Q
            index_M=svm_struct.index_M;
            local_Mx=x(index_M,:);
            local_qM=SemiSVM.Kernel(x,local_Mx);
            local_qQ=SemiSVM.Kernel(x,x(index_Q,:));
            gama=local_qM*beta(2:end,:)+beta(1,:)+local_qQ*y(index_Q)*d;
        end
        function    beta=SPGetBeta(x,y,svm_struct,Qc)
            Q=IncSemiSVM_algorithm2.SPFindQR(x,svm_struct);
            
            if det(Q)<1e-42
                m=length(Q);
                n=eye(m);
                Q=Q+1e-7*n;
            end
            beta=linsolve(Q,Qc);
            
            %beta=mldivide(Q,Qc);
            %beta=linsolve(Q,Qc);
        end
        function Q=SPFindQR(x,svm_struct)
            index_M=svm_struct.index_M;
            length_M=length(index_M);
            x_M=x(index_M,:);
            Q=SemiSVM.Kernel(x_M,x_M);
            tmp1=ones(length_M,1);
            tmp1p=[0;tmp1];
            Q=[tmp1';Q];
            Q=[tmp1p  Q];
        end
        function   [out_g,out_svm_struct,out_Qc,direction]=SPUpdate...
                (svm_struct,g,beta,gama,max_delta,index,x,y,direction,index_class)
            global index_Q flag C
            [out_svm_struct,direction,index_A]=IncSemiSVM_algorithm2.SPUpdataSvmstruct(svm_struct,beta,max_delta,index,y,direction,index_class);
            [out_g,pms]=IncSemiSVM_algorithm2.SPUpdataG(g,gama,max_delta,svm_struct);
            % other case
            if ~isempty(index_A)
                fprintf('A ~isempty \n');
                
                %[out_svm_struct,out_g]=CCCP.Run(x,y,out_svm_struct,out_g,index_A);
                flag=1;
            end
            out_svm_struct.pms=pms;
            index_Q=[out_svm_struct.index_E(out_svm_struct.alpha(out_svm_struct.index_E)~=0 & out_svm_struct.alpha(out_svm_struct.index_E)~=C);out_svm_struct.index_O(out_svm_struct.alpha(out_svm_struct.index_O)~=0)];
            [out_Qc]=IncSemiSVM_algorithm2.SPGetQc(out_svm_struct,direction);
        end
        function [out_g,pms]=SPUpdataG(g,gama,max_delta,out_svm_struct)
            global C C_star
            other_set=[1:length(g)];
            other_set(out_svm_struct.index_M)=[];
            g(other_set)=g(other_set)+gama(other_set)*max_delta;
         
            y=out_svm_struct.y;
            label_size=out_svm_struct.label_size;
            LowBoud=0.*(y(1:label_size)>0)+(-C).*(y(1:label_size)<0);
            UppBound=(C).*(y(1:label_size)>0)+0.*(y(1:label_size)<0);
           % gi=out_svm_struct.gi;
            y([1:label_size])=[];
            gi=g;
            gi([1:label_size])=[];
            os=out_svm_struct;
            local_f=(os.b-gi).*y+1;
            mu=(C-C_star)*(local_f<0);
            low_tmp=(-mu).*(y>0)+(mu-C+C_star).*(y<0);
            upp_tmp=(C-C_star-mu).*(y>0)+(mu).*(y<0);
            LowBound=[LowBoud;low_tmp];
            UppBound=[UppBound;upp_tmp];
            pms=[LowBound UppBound];
            out_g=g;
        end
        
        function  [out_svm_struct,direction_A,index_A] = SPUpdataSvmstruct...
                (svm_struct,beta,max_delta,index,y,d,index_class)
            global  fake_zero C C_end C_star flag
            index_A=[];
            local_beta=beta;
            local_beta_length=length(local_beta);
            beta_b=beta(1);
            beta_s=beta(2:local_beta_length);    
            alpha=svm_struct.alpha;
            support_set=svm_struct.index_M;
            b=svm_struct.b;  
            alpha(support_set)=alpha(support_set)+beta_s*max_delta;
            b=b+beta_b*max_delta;
            svm_struct.b=b;
%             k1=(svm_struct.pms(:,2)>C-fake_zero);
%             k2=(svm_struct.pms(:,1)<-C+fake_zero);
            C_star=C_star+max_delta*d;
%             svm_struct.pms(k1,2)=C;
%             svm_struct.pms(k2,1)=-C;
%             alpha(svm_struct.index_E)=svm_struct.pms(svm_struct.index_E,2);
%             alpha(svm_struct.index_O)=svm_struct.pms(svm_struct.index_O,1);
%             svm_struct.alpha=alpha;
            direction_A=C_end-C_star;
            fprintf('delta=%.4d,C_star=%.4d\n',max_delta*d,C_star);
            switch  index_class
                case 1
                    tmp_index = svm_struct.index_M == index;
                    svm_struct.index_M(tmp_index)=[];
                    svm_struct.index_O = sort([svm_struct.index_O;index]);
                case 2
                    tmp_index = svm_struct.index_M == index;
                    svm_struct.index_M(tmp_index)=[];
                    svm_struct.index_E = sort([svm_struct.index_E;index]);
                case 3
                    tmp_index = svm_struct.index_O == index;   
                    svm_struct.index_O(tmp_index)=[];
                    svm_struct.index_M = sort([svm_struct.index_M;index]);
                case 4
                    tmp_index = svm_struct.index_E == index;
                    svm_struct.index_E(tmp_index)=[];
                    svm_struct.index_M = sort([svm_struct.index_M;index]);
                case 8
                    tmp_index = svm_struct.index_O == index;
                    svm_struct.index_O(tmp_index)=[];
                    index_A = index;
                    if svm_struct.pms(index,1)<0
                        svm_struct.pms(index,1)=0;
                        svm_struct.pms(index,2)=C-C_star;
                    else
                        svm_struct.pms(index,1)=-C+C_star;
                        svm_struct.pms(index,2)=0;                       
                    end
                case 9
                    tmp_index = svm_struct.index_E == index;
                    svm_struct.index_E(tmp_index)=[];
                    index_A = index;
                    if svm_struct.pms(index,1)<0
                        svm_struct.pms(index,1)=0;
                        svm_struct.pms(index,2)=C-C_star;
                    else
                        svm_struct.pms(index,1)=-C+C_star;
                        svm_struct.pms(index,2)=0;                     
                    end
                otherwise
                    flag=1;
            end
            out_svm_struct=svm_struct;
        end
 
        function  [max_delta_zeta,index,index_class]=SPGetMaxChange...
                (svm_struct,g,beta,gama,y,direction,original_y)
            [max_change_alpha,index_max_change_alpha,class_alpha]=IncSemiSVM_algorithm2.SPFindChangeAlpha(svm_struct,beta,direction);
            [max_change_g,index_max_change_g,class_g]=IncSemiSVM_algorithm2.SPFindMaxChangeG(svm_struct,g,gama);
            [max_change_toA,index_max_change_toA,class_toA]=IncSemiSVM_algorithm2.SPFindChangeToA...
                (svm_struct,y,g,gama,original_y);
            [max_delta_zeta,index,index_class]=IncSemiSVM_algorithm2.SPFindMaxChange...
                (max_change_alpha, index_max_change_alpha,max_change_g,...
                index_max_change_g,class_alpha,class_g,max_change_toA,index_max_change_toA,class_toA);
        end
        function    [updata_max_change, updata_index,index_class]=SPFindMaxChange...
                (max_change_alpha, index_max_change_alpha,max_change_g,...
                index_max_change_g,class_alpha,class_g,max_change_toA,index_max_change_toA,class_toA)
            updata_max_change=inf;
            updata_index=0;
            index_class=inf;
            if max_change_alpha <= updata_max_change
                updata_max_change=max_change_alpha;
                updata_index=index_max_change_alpha;
                index_class=class_alpha;
            end
            if max_change_g  <= updata_max_change 
                updata_max_change=max_change_g;
                updata_index=index_max_change_g;
                index_class=class_g;
            end
            if max_change_toA <=updata_max_change
                updata_max_change=max_change_toA;
                updata_index=index_max_change_toA;
                index_class=class_toA;
            end
        end

        function [max_change_toA,index_max_change_toA, class_toA]=SPFindChangeToA...
                (svm_struct,y,g,gama,original_y)
            global fake_zero
            class_toA = inf;
            other_set=[1:length(y)];
            local_g=g;
            index_max_change_toA=0;
            local_index=find(gama<fake_zero & gama>-fake_zero);
            gama(local_index)=[];
            other_set(local_index)=[];
            local_g(local_index)=[];
            local_g=local_g.*y(other_set)+1;
            gama=gama.*y(other_set);            
            local_index3=find(local_g.*gama>=0);
            local_g(local_index3)=[];
            gama(local_index3)=[];
            other_set(local_index3)=[];
            local_change_g = -local_g./gama;
            [local_change, local_index]=min(local_change_g);
            if isempty(local_change_g)
                max_change_toA=inf;
            else
                max_change_toA=local_change;
                index_max_change_toA=other_set(local_index);
                    if g(index_max_change_toA) < 0
                        class_toA = 9;  % g moves from E set to A set
                    else
                        class_toA = 8;   % sample moves from O set to A set
                    end
            end
        end
        
        function [max_change_alpha,index_max_change_alpha,class_alpha] = SPFindChangeAlpha(svm_struct,beta,d)
            global fake_zero 
            index_max_change_alpha=0;
            class_alpha=inf;
            support_set=svm_struct.index_M;
            alpha=svm_struct.alpha;
            local_y=svm_struct.y(support_set);
            local_alpha=alpha(support_set);
            local_beta_length=length(beta);
            local_beta=beta(2:local_beta_length);
            local_index=find(local_beta<fake_zero & local_beta>-fake_zero);
            local_alpha(local_index)=[];
            support_set(local_index)=[];
            local_beta(local_index)=[];
            local_y(local_index)=[];
            local_Cdown=svm_struct.pms(support_set,1).*(local_beta<0);
            local_change_aldown=(local_Cdown-local_alpha)./local_beta;
            local_Cup=svm_struct.pms(support_set,2).*(local_beta>0);
            local_change_alup=(local_Cup-local_alpha)./(local_beta-d);
            local_change_al=local_change_aldown.*(local_beta<0)+local_change_alup.*(local_beta>0);
            local_Cdown1=svm_struct.pms(support_set,1).*(local_beta<0);
            local_change_aldown1=(local_Cdown1-local_alpha)./(local_beta+d);
            local_Cup1=svm_struct.pms(support_set,2).*(local_beta>0);
            local_change_alup1=(local_Cup1-local_alpha)./local_beta;
            local_change_al1=local_change_aldown1.*(local_beta<0)+local_change_alup1.*(local_beta>0);
            local_change_al(local_change_al<=0)=inf;
            local_change_al1(local_change_al1<=0)=inf;
            local_change_al2=local_change_al.*(local_y>0)+local_change_al1.*(local_y<0);
            [min_alpha_change,min_alpha_change_index]=min(local_change_al2);
            if isempty(local_change_al2)
                max_change_alpha=inf;
            else
                max_change_alpha=min_alpha_change;
                index_max_change_alpha=support_set(min_alpha_change_index);
                if local_beta(min_alpha_change_index)>0
                    class_alpha=2;  % alpha moves from M set to E set
                else
                    class_alpha=1;  % alpha moves from M set to O set
                end
            end
        end
        function [max_change_g,index_max_change_g,class_g] = SPFindMaxChangeG(svm_struct,g,gama)
            class_g = inf;
            other_set=[1:length(g)];
            other_set(svm_struct.index_M)=[];
            gama=gama(other_set);       
            local_g=g(other_set);
            index_max_change_g=0;
            local_index3=find(local_g.*gama>=0);
            local_g(local_index3)=[];
            gama(local_index3)=[];
            other_set(local_index3)=[];
            local_change_g=-local_g./gama;
             local_change_g(local_change_g<1e-10)=inf;
            [local_change,local_index]=min(local_change_g);
            if isempty(local_change_g)
                max_change_g=inf;
            else
                max_change_g=local_change;
                index_max_change_g=other_set(local_index);
                if gama(local_index)>0
                    class_g=4;  % g moves from E set to M set
                else
                    class_g=3;  % g moves from O set to M set
                end
            end
        end
    end    
end