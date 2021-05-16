classdef CCCP
    properties
        KKT=0;     
    end
    
    methods(Static = true)
        function obj=IncSemiSVM(a)
            obj.KKT=a;
        end
        
        function [out,g]=Run(original_x,original_y,out_svm_struct,out_g,index)
         
            g=out_g;
             svm_struct=out_svm_struct;
                [svm_struct,index_A,direction_A]=CCCP.InitialA(svm_struct,g(index),index);
                [Qc]=CCCP.GetQc(svm_struct,original_x,index_A,direction_A);
                while ~isempty(index_A)
                    if direction_A==0
                        break;
                    end
                    [beta]=CCCP.GetBeta(original_x,original_y,svm_struct,Qc);
                    [gama] = CCCP.GetGama(original_x,original_y,svm_struct,beta,direction_A,index_A);
                    [max_delta,index,index_class]=CCCP.GetMaxChange...
                        (svm_struct,g,beta,gama,original_y,direction_A,index_A,original_y);
                    [g,svm_struct,Qc,direction_A,index_A]=CCCP.Update...
                        (svm_struct,g,beta,gama,max_delta,index,original_x,original_y,direction_A,index_A,index_class);
                %    steps=steps+1;
             
                end
       
        
            out=svm_struct;
                M=g(svm_struct.index_M);
              E=g(svm_struct.index_E);
              O=g(svm_struct.index_O);
        end
        
        function   [out_g,out_svm_struct,out_Qc,direction_A,index_A]=Update...
                (svm_struct,g,beta,gama,max_delta,index,x,y,direction_A,index_A,index_class)
            [out_svm_struct,direction_A,index_A]=CCCP.UpdataSvmstruct(svm_struct,beta,max_delta,index,y,direction_A,index_A,index_class);
            [out_g]=CCCP.UpdataG(g,svm_struct,gama,max_delta);
            [out_Qc]=CCCP.GetQc(out_svm_struct,x,index_A,direction_A);
        end
        
        function [out_g]=UpdataG(g,svm_struct,gama,max_delta)
            other_set=[1:length(g)];
            other_set(svm_struct.index_M)=[];
            g(other_set)=g(other_set)+gama(other_set)*max_delta;
            out_g=g;
        end
        function  [out_svm_struct,direction_A,index_A] = UpdataSvmstruct...
                (svm_struct,beta,max_delta,index,y,direction_A,index_A,index_class)
            global C C_star
            local_beta=beta;
            local_beta_length=length(local_beta);
            beta_b=beta(1);
            beta_s=beta(2:local_beta_length);
            alpha=svm_struct.alpha;
            support_set=svm_struct.index_M;
            b=svm_struct.b;
            alpha(support_set)=alpha(support_set)+beta_s*max_delta;
            alpha(index_A) = alpha(index_A)+direction_A*max_delta;
        %    direction_A=direction_A-direction_A*max_delta;
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
                    tmp_index = svm_struct.index_M == index;
%                        if tmp_index==0
%                         k=1;
%                     end
                    svm_struct.index_M(tmp_index)=[];
                    svm_struct.index_O = sort([svm_struct.index_O;index]);
                case 2
                    tmp_index = svm_struct.index_M == index;
%                        if tmp_index==0
%                         k=1;
%                     end
                    svm_struct.index_M(tmp_index)=[];
                    svm_struct.index_E = sort([svm_struct.index_E;index]);
                case 3
                    tmp_index = svm_struct.index_O == index;
%                        if tmp_index==0
%                         k=1;
%                     end
                    svm_struct.index_O(tmp_index)=[];
                    svm_struct.index_M = sort([svm_struct.index_M;index]);
                case 4
                    tmp_index = svm_struct.index_E == index;
%                        if tmp_index==0
%                         k=1;
%                     end
                    svm_struct.index_E(tmp_index)=[];
                    svm_struct.index_M = sort([svm_struct.index_M;index]);
                case 5
                    tmp_index = find(index_A == index);
%                        if tmp_index==0
%                         k=1;
%                     end
                    index_A(tmp_index)=[];
                    direction_A(tmp_index,:) = [];
                    svm_struct.index_M = sort([svm_struct.index_M;index]);
                case 6
                    tmp_index = find(index_A == index);
%                        if tmp_index==0
%                         k=1;
%                     end
                    index_A(tmp_index)=[];
                    direction_A(tmp_index,:) = [];
                    svm_struct.index_O = sort([svm_struct.index_O;index]);
                case 7
                    tmp_index = find(index_A == index);
%                        if tmp_index==0
%                         k=1;
%                     end
                    index_A(tmp_index)=[];
                    direction_A(tmp_index,:) = [];
                    svm_struct.index_E = sort([svm_struct.index_E;index]);
                case 8
                    tmp_index = find(svm_struct.index_O == index);
%                        if tmp_index==0
%                         k=1;
%                     end
                    svm_struct.index_O(tmp_index)=[];
                    index_A = [index_A;index];
                    if svm_struct.pms(index,1)<0
                        svm_struct.pms(index,1)=0;
                        svm_struct.pms(index,2)=C-C_star;
                        tmp1=C-C_star;
                    else
                        svm_struct.pms(index,1)=-C+C_star;
                        svm_struct.pms(index,2)=0;
                        tmp1=-C+C_star;                        
                    end
                    local_direction = tmp1-svm_struct.alpha(index);
                    direction_A = [direction_A;local_direction];
                case 9
                    tmp_index = find(svm_struct.index_E == index);
%                        if tmp_index==0
%                         k=1;
%                     end
                    svm_struct.index_E(tmp_index)=[];
                    index_A = [index_A; index];
                    if svm_struct.pms(index,1)<0
                        svm_struct.pms(index,1)=0;
                        svm_struct.pms(index,2)=C-C_star;
                        tmp1=C-C_star;
                    else
                        svm_struct.pms(index,1)=-C+C_star;
                        svm_struct.pms(index,2)=0;
                        tmp1=-C+C_star;                        
                    end
                    local_direction = tmp1-svm_struct.alpha(index);
                    direction_A = [direction_A;local_direction];
                case 10
                    if svm_struct.pms(index,1)<0
                        svm_struct.pms(index,1)=0;
                        svm_struct.pms(index,2)=C-C_star;
                        tmp1=C-C_star;
                    else
                        svm_struct.pms(index,1)=-C+C_star;
                        svm_struct.pms(index,2)=0;
                        tmp1=-C+C_star;                        
                    end
                    local_direction = tmp1-svm_struct.alpha(index);
                    tmp2=find(index_A==index);
                    direction_A(tmp2) = local_direction;
           
                   
            end
            out_svm_struct=svm_struct;
        end
 
        function  [max_delta_zeta,index,index_class]=GetMaxChange...
                (svm_struct,g,beta,gama,y,direction_A,index_A,original_y)
            [max_change_alpha,index_max_change_alpha,class_alpha]=CCCP.FindChangeAlpha(svm_struct,beta);
            [max_change_g,index_max_change_g,class_g]=CCCP.FindMaxChangeG(svm_struct,g,gama);
            [max_change_A,index_max_change_A,class_A]=CCCP.FindChangeA...
                (svm_struct,y,beta, direction_A,index_A,g,gama);
            [max_change_toA,index_max_change_toA,class_toA]=CCCP.FindChangeToA...
                (svm_struct,y, index_A,g,gama,original_y);
            [max_delta_zeta,index,index_class]=CCCP.FindMaxChange...
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
            [min_alpha_change  min_alpha_change_index]=min(local_change_al);
            if min_alpha_change < max_change_A
                max_change_A = min_alpha_change;
                index_max_change_A=local_index_A(min_alpha_change_index);
                class_A=6;  % alpha moves from A set to O set
            end             
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
        function [max_change_toA,index_max_change_toA, class_toA]=FindChangeToA...
                (svm_struct,y, index_A,g,gama,original_y)
            % condition 5 in the paper
            global fake_zero
            class_toA = inf;
            other_set=[1:length(y)];
            local_g=g;
            gama=gama;
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
                 local_change_g(local_change_g<100*fake_zero)=inf;
            [local_change, local_index]=min(local_change_g);
            if length(local_change_g)==0
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
        
        function [max_change_alpha,index_max_change_alpha,class_alpha] = FindChangeAlpha(svm_struct,beta)
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
            local_change_al(local_change_al<0)=inf;
            [min_alpha_change,min_alpha_change_index]=min(local_change_al);
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
            gama=gama(other_set);        
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
             local_change_g(local_change_g<1e-8)=inf;
            [local_change,local_index]=min(local_change_g);
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
        function [svm_struct,index_A,direction_A] = InitialA(svm_struct,gc,index)
            
            index_A=[index];
            pms=svm_struct.pms(index,:);
            y=svm_struct.y(index);
            if svm_struct.alpha(index)==0
                direction_A=(y>0)*pms(:,2)+(y<0)*pms(:,1);
            else                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                direction_A=(y<0)*pms(:,2)+(y>0)*pms(:,1);
            end
%             svm_struct.alpha(index)=[];
%             svm_struct.pms(index,:)=[];
%             svm_struct.y(index)=[];
%             svm_struct.x(index)=[];
        end
        
        function [Qc]=GetQc(svm_struct,x,index_A,direction_A)
            Qc=[];
            local_direction=direction_A;
            if ~isempty(index_A)
                index_M=svm_struct.index_M;
                length_A=length(index_A);
                tmp1=ones(1,length_A);
                x_M=x(index_M,:);
                x_A=x(index_A,:);
                H_MA=SemiSVM.Kernel(x_M,x_A);
                tmp_H=[tmp1;H_MA];
                Qc=-tmp_H*local_direction;
            end
        end
        function [gama]=GetGama(x,y,svm_struct,beta,direction_A,index_A)
            support_set=svm_struct.index_M;
            local_x=x;
            local_sx=x(support_set,:);
            local_Ax=x(index_A,:);
            local_qbeta=SemiSVM.Kernel(local_x,local_sx);
            local_qA=SemiSVM.Kernel(local_x,local_Ax);
            local_qbeta=[ones(length(y),1)  local_qbeta];
            gama=local_qbeta*beta+local_qA*direction_A;
        end
        function beta=GetBeta(x,y,svm_struct,Qc)
            Q=CCCP.FindQR(x,y,svm_struct);
           if det(Q)<1e-42
                m=length(Q);
                n=eye(m);
                Q=Q+1e-7*n;
            end
            beta=linsolve(Q,Qc);
        end
        function Q=FindQR(x,y,svm_struct)
            index_M=svm_struct.index_M;
            length_M=length(index_M);
            x_M=x(index_M,:);
            Q=SemiSVM.Kernel(x_M,x_M);
            tmp1=ones(length_M,1);
            tmp1p=[0;tmp1];
            Q=[tmp1';Q];
            Q=[tmp1p  Q];
        end
        
    end
    
end

