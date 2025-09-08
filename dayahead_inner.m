function dayahead_inner(Y_syslim)
%% 形参说明
% Y_syslim为系统总超调抑制能力下限
%% 定义辅助常数
my_Phi = 0.5;
%% 数据及参数导入
load dataimport
%% 火电机组变量
u_G = binvar(dhdata.T,Gpara.N,'full');          %机组启停状态，开机为1，关机为0
P_G = sdpvar(dhdata.T,Gpara.N,'full');          %机组出力
A_G = sdpvar(dhdata.T,Gpara.N,Gpara.N,'full');  %辅助变量，A_G(t,m,n) = u_G(t,m)*P_G(t,n)
%% 电化学储能变量
u_ES = binvar(dhdata.T,ESpara.N,'full');        %储能工作状态，放电为1，充电为0
P_ES = sdpvar(dhdata.T,ESpara.N,'full');        %储能输出功率，放电为正，充电为负
E_ES = sdpvar(dhdata.T,ESpara.N,'full');        %储能剩余能量
A_ES = sdpvar(dhdata.T,Gpara.N,'full');         %辅助变量，A_ES(t,m) = u_G(t,m)*P_ES(t)
%% 风电变量
P_wind = sdpvar(dhdata.T,N_wind,'full');        %风电实际出力
%% 目标函数
C_G_run = 0.25*sum(sum(repmat(Gpara.a,dhdata.T,1).*P_G+repmat(Gpara.b,dhdata.T,1).*u_G)); %火电机组运行成本
C_G_on = sum(sum(repmat(Gpara.on,dhdata.T-1,1).*(u_G(2:dhdata.T,:)-u_G(1:dhdata.T-1,:)+abs(u_G(2:dhdata.T,:)-u_G(1:dhdata.T-1,:)))/2)); %火电机组启停成本
C_wind = 0.25*sum(sum(penalty_wind*(dhdata.wind-P_wind)));   %弃风成本
obj = C_G_run+C_G_on+C_wind;
%% 火电机组约束条件
constraints = [];
% 火电机组出力约束
constraints = [constraints,u_G.*repmat(Gpara.Pmin,dhdata.T,1)<=P_G<=u_G.*repmat(Gpara.Pmax,dhdata.T,1)];
% 火电机组爬坡率约束
constraints = [constraints,-repmat(15*Gpara.ramp,dhdata.T-1,1)-(2-u_G(2:dhdata.T,:)-u_G(1:dhdata.T-1,:))*M<=...
    P_G(2:dhdata.T,:)-P_G(1:dhdata.T-1,:)<=repmat(15*Gpara.ramp,dhdata.T-1,1)+(2-u_G(2:dhdata.T,:)-u_G(1:dhdata.T-1,:))*M];
% 火电机组启停时间约束
for n = 1:Gpara.N
    for t = 1:dhdata.T-4*Gpara.sst(n)
        constraints = [constraints,sum(u_G(t+1:t+4*Gpara.sst(n),n))>=4*Gpara.sst(n)*(u_G(t+1,n)-u_G(t,n))];         %最小开机时间约束(开关机时间外)
        constraints = [constraints,sum(1-u_G(t+1:t+4*Gpara.sst(n),n))>=4*Gpara.sst(n)*(u_G(t,n)-u_G(t+1,n))];       %最小停机时间约束(开关机时间外)
    end
    for t = dhdata.T-4*Gpara.sst(n)+1:dhdata.T-1
        constraints = [constraints,sum(u_G(t+1:dhdata.T,n))>=(dhdata.T-t)*(u_G(t+1,n)-u_G(t,n))];                   %最小开机时间约束(开关机时间内)
        constraints = [constraints,sum(1-u_G(t+1:dhdata.T,n))>=(dhdata.T-t)*(u_G(t,n)-u_G(t+1,n))];                 %最小停机时间约束(开关机时间内)
    end
end
% % 常规波动下火电机组一次调频稳态功率变化量约束 P^{G,L}_{t,i}
% for n = 1:Gpara.N
%     for t = 1:dhdata.T
%         constraints = [constraints,Gpara.Pmax(n)*(dhdata.el(t)+sum(dhdata.ew(t,:)))/Gpara.adj(n)+...
%             (Gpara.Pmax(n)*df.deadzone/Gpara.adj(n)-Gpara.Pmax(n))*(k_load*dhdata.load(t)+ESpara.Pmax*ESpara.kd+sum(u_G(t,:).*Gpara.Pmax./Gpara.adj))<=...
%             -P_G(t,n)*(k_load*dhdata.load(t)+ESpara.Pmax*ESpara.kd)-sum(A_G(t,:,n).*Gpara.Pmax./Gpara.adj)+(1-u_G(t,n))*M];
%         constraints = [constraints,Gpara.Pmax(n)*(dhdata.el(t)+sum(dhdata.ew(t,:)))/Gpara.adj(n)+...
%             (Gpara.Pmax(n)*df.deadzone/Gpara.adj(n)+Gpara.Pmin(n))*(k_load*dhdata.load(t)+ESpara.Pmax*ESpara.kd+sum(u_G(t,:).*Gpara.Pmax./Gpara.adj))<=...
%             P_G(t,n)*(k_load*dhdata.load(t)+ESpara.Pmax*ESpara.kd)+sum(A_G(t,:,n).*Gpara.Pmax./Gpara.adj)+(1-u_G(t,n))*M];
%         constraints = [constraints,repmat(P_G(t,n),1,Gpara.N)-(1-u_G(t,:))*M<=A_G(t,:,n)<=repmat(P_G(t,n),1,Gpara.N)+(1-u_G(t,:))*M];
%         constraints = [constraints,-u_G(t,:)*M<=A_G(t,:,n)<=u_G(t,:)*M];
%     end
% end
% new 火电一次调频约束(LQ1 paper)
for n = 1: Gpara.N
    for t = 1:dhdata.T
        constraints = [constraints,Gpara.kp(n).*sum(A_G(t,:,n).*Gpara.Pmax./Gpara.adj) + P_G(t, n).*(Gpara.kp(n)*k_load*dhdata.load(t) + Gpara.kp(n)*ESpara.Pmax*ESpara.kd-Gpara.kp(n))-...
            (dhdata.el(t)+sum(dhdata.ew(t,:)+my_Phi).*(1+Gpara.adj(n).*Gpara.kp(n))) <= ...
            Gpara.Pmax(n).*(Gpara.kp(n).*sum(u_G(t,:).*Gpara.Pmax./Gpara.adj)+(Gpara.Pmax(n).*k_load*dhdata.load(t)+Gpara.kp(n).*ESpara.Pmax.*ESpara.kd-Gpara.kp)) + (1-u_G(t,n))*M];
        constraints = [constraints,Gpara.kp(n)*sum(A_G(t,:,n).*Gpara.Pmax./Gpara.adj) + P_G(t, n).*(Gpara.kp(n).*k_load*dhdata.load(t) + Gpara.kp(n).*ESpara.Pmax.*ESpara.kd-Gpara.kp(n))+...
            (dhdata.el(t)+sum(dhdata.ew(t,:)+my_Phi).*(1+Gpara.adj(n)*Gpara.kp(n))) >= ...
            Gpara.Pmin(n).*(Gpara.kp(n).*sum(u_G(t,:).*Gpara.Pmax./Gpara.adj)+(Gpara.Pmax(n).*k_load.*dhdata.load(t)+Gpara.kp(n).*ESpara.Pmax.*ESpara.kd-Gpara.kp)) + (1-u_G(t,n))*M];
         constraints = [constraints,repmat(P_G(t,n),1,Gpara.N)-(1-u_G(t,:))*M<=A_G(t,:,n)<=repmat(P_G(t,n),1,Gpara.N)+(1-u_G(t,:))*M];
         constraints = [constraints,-u_G(t,:)*M<=A_G(t,:,n)<=u_G(t,:)*M];
    end
end
%% 风电机组约束条件
constraints = [constraints,0<=P_wind<=dhdata.wind];
%% 电化学储能约束条件
% 电化学储能功率约束
constraints = [constraints,-(1-u_ES)*M<=P_ES<=ESpara.Pmax+(1-u_ES)*M];
constraints = [constraints,-ESpara.Pmax-u_ES*M<=P_ES<=u_ES*M];
% 电化学储能能量约束
constraints = [constraints,E_ES(1)==ESpara.Eini];
constraints = [constraints,ESpara.Emin<=E_ES<=ESpara.Emax];
constraints = [constraints,0.25*P_ES(1:dhdata.T-1)/ESpara.eff-(1-u_ES(1:dhdata.T-1))*M<=E_ES(1:dhdata.T-1)-E_ES(2:dhdata.T)<=0.25*P_ES(1:dhdata.T-1)/ESpara.eff+(1-u_ES(1:dhdata.T-1))*M];
constraints = [constraints,0.25*P_ES(1:dhdata.T-1)*ESpara.eff-u_ES(1:dhdata.T-1)*M<=E_ES(1:dhdata.T-1)-E_ES(2:dhdata.T)<=0.25*P_ES(1:dhdata.T-1)*ESpara.eff+u_ES(1:dhdata.T-1)*M];
constraints = [constraints,ESpara.Eini-(1-u_ES(dhdata.T))*M<=E_ES(dhdata.T)-0.25*P_ES(dhdata.T)/ESpara.eff<=ESpara.Emax+(1-u_ES(dhdata.T))*M];
constraints = [constraints,ESpara.Eini-u_ES(dhdata.T)*M<=E_ES(dhdata.T)-0.25*P_ES(dhdata.T)*ESpara.eff<=ESpara.Emax+u_ES(dhdata.T)*M];
% % 常规波动下电化学储能一次调频稳态功率变化量约束
% for t = 1:dhdata.T
%     constraints = [constraints,ESpara.Pmax*ESpara.kd*(dhdata.el(t)+sum(dhdata.ew(t,:)))+...
%         ESpara.Pmax*(ESpara.kd*df.deadzone-1)*(k_load*dhdata.load(t)+ESpara.Pmax*ESpara.kd+sum(u_G(t,:).*Gpara.Pmax./Gpara.adj))<=...
%         -P_ES(t)*(k_load*dhdata.load(t)+ESpara.Pmax*ESpara.kd)-sum(A_ES(t,:).*Gpara.Pmax./Gpara.adj)];
%     constraints = [constraints,ESpara.Pmax*ESpara.kd*(dhdata.el(t)+sum(dhdata.ew(t,:)))+...
%         ESpara.Pmax*(ESpara.kd*df.deadzone-1)*(k_load*dhdata.load(t)+ESpara.Pmax*ESpara.kd+sum(u_G(t,:).*Gpara.Pmax./Gpara.adj))<=...
%         P_ES(t)*(k_load*dhdata.load(t)+ESpara.Pmax*ESpara .kd)+sum(A_ES(t,:).*Gpara.Pmax./Gpara.adj)];
%     constraints = [constraints,repmat(P_ES(t),1,Gpara.N)-(1-u_G(t,:))*M<=A_ES(t,:)<=repmat(P_ES(t),1,Gpara.N)+(1-u_G(t,:))*M];
%     constraints = [constraints,-u_G(t,:)*M<=A_ES(t,:)<=u_G(t,:)*M];
% end
% new 电化学储能一次调频约束(LQ1 paper)类似上述火电
for t = 1:dhdata.T
    constraints = [constraints, sum(A_ES(t,:).*Gpara.Pmax./Gpara.adj) + P_ES(t)*(k_load*dhdata.load(t) + ESpara.Pmax*ESpara.kd-1) -...
        (dhdata.el(t)+sum(dhdata.ew(t,:)+my_Phi)*ESpara.kv) <=...
        ESpara.Pmax*(sum(u_G(t,:).*Gpara.Pmax./Gpara.adj) + k_load*dhdata.load(t) + ESpara.Pmax*ESpara.kd-1)];
    constraints = [constraints, sum(A_ES(t,:).*Gpara.Pmax./Gpara.adj) + P_ES(t)*(k_load*dhdata.load(t) + ESpara.Pmax*ESpara.kd-1) +...
        (dhdata.el(t)+sum(dhdata.ew(t,:)+my_Phi)*ESpara.kv) >=...
        -ESpara.Pmax*(sum(u_G(t,:).*Gpara.Pmax./Gpara.adj) + k_load*dhdata.load(t) + ESpara.Pmax*ESpara.kd-1)];
    constraints = [constraints,repmat(P_ES(t),1,Gpara.N)-(1-u_G(t,:))*M<=A_ES(t,:)<=repmat(P_ES(t),1,Gpara.N)+(1-u_G(t,:))*M];
    constraints = [constraints,-u_G(t,:)*M<=A_ES(t,:)<=u_G(t,:)*M];
end
%% 系统运行约束条件
% 系统功率平衡约束
constraints = [constraints,sum(P_G,2)+sum(P_wind,2)+P_ES==dhdata.load];
% 常规波动下的稳态频差约束
constraints = [constraints,(dhdata.el+sum(dhdata.ew,2))/(df.steady1-df.deadzone)-k_load*dhdata.load<=ESpara.Pmax*ESpara.kd+sum(u_G.*repmat(Gpara.Pmax./Gpara.adj,dhdata.T,1),2)];
% 常规波动下的超调频差约束(超调抑制能力约束)
constraints = [constraints,sum(u_G.*repmat(Gpara.Y,dhdata.T,1),2)>=Y_syslim];
% 常规波动下的频率变化率约束
constraints = [constraints,(dhdata.el+sum(dhdata.ew,2))/df.changerate1<=ESpara.Pmax*ESpara.kv+2*sum(u_G.*repmat(Gpara.Pmax.*Gpara.H,dhdata.T,1),2)];
% 系统正备用约束
for n = 1:Gpara.N
    constraints = [constraints,sum(u_G.*repmat(Gpara.Pmax,dhdata.T,1)-P_G,2)-(u_G(:,n)*Gpara.Pmax(n)-P_G(:,n))+(ESpara.Pmax-P_ES)>=dhdata.el+sum(dhdata.ew,2)+P_G(:,n)];
    constraints = [constraints,sum(u_G.*repmat(Gpara.Pmax,dhdata.T,1)-P_G,2)-(u_G(:,n)*Gpara.Pmax(n)-P_G(:,n))+(4*ESpara.eff*(E_ES-ESpara.Emin)-P_ES)>=dhdata.el+sum(dhdata.ew,2)+P_G(:,n)];
end
% 系统负备用约束
constraints = [constraints,sum(P_G-u_G.*repmat(Gpara.Pmin,dhdata.T,1),2)+(P_ES+ESpara.Pmax)>=dhdata.el+sum(dhdata.ew,2)];
constraints = [constraints,sum(P_G-u_G.*repmat(Gpara.Pmin,dhdata.T,1),2)+(P_ES-4*(E_ES-ESpara.Emax)/ESpara.eff)>=dhdata.el+sum(dhdata.ew,2)];
% new 跌落最低点约束(LQ1 paper)


%% 求解
ops = sdpsettings('solver','cplex','verbose',1);
% ops.gurobi.TimeLimit = 200;
% ops.gurobi.MIPGap = 0.001;
optimize(constraints,obj,ops)
%% 结果保存
% 火电机组优化结果
u_G_dh = value(u_G);
P_G_dh = value(P_G);
% 风电机组优化结果
P_wind_dh = value(P_wind);
% 电化学储能优化结果
u_ES_dh = value(u_ES);
P_ES_dh = value(P_ES);
E_ES_dh = value(E_ES);
% 成本优化结果
C_G_run_dh = value(C_G_run);
C_G_on_dh = value(C_G_on);
C_wind_dh = value(C_wind);
C_sum_dh = value(obj);
% 保存
save dh_result u_G_dh P_G_dh P_wind_dh u_ES_dh P_ES_dh E_ES_dh C_G_run_dh C_G_on_dh C_wind_dh C_sum_dh
end