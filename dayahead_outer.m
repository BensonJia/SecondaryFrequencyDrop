clear;clc;tic
%% 数据及参数导入
DataImport
Y_syslim = zeros(dhdata.T,1);   %系统超调抑制能力下限初值
times = 0;                      %迭代次数
my_Phi = 0.5;                   %辅助常量
t_1 = 14;
%% 迭代
while 1
    str=['Iterating times: ' num2str(times)];
    disp(str);
    % 日前调度
    dayahead_inner(Y_syslim);
    times = times+1;
    load dh_result
    % 火电机组等效参数求解
    adj_sys = 1./sum(u_G_dh.*repmat(Gpara.Pmax./Gpara.adj,dhdata.T,1),2);   %等效调差系数
    innerpara1 = sum(u_G_dh.*repmat(Gpara.Pmax./Gpara.adj.*(1-1./(1+Gpara.adj.*Gpara.kp)),dhdata.T,1),2);
    kp_sys = innerpara1./(1-adj_sys.*innerpara1);   %等效比例系数
    t_st = -(1+Gpara.adj.*Gpara.kp).*log(0.02*(1+(1+Gpara.adj.*Gpara.kp)))./(Gpara.adj.*Gpara.ki);
    t_mid = 0.5*adj_sys.*sum(u_G_dh.*Gpara.Pmax./Gpara.adj.*t_st,2);
    innerpara2 = sum(u_G_dh.*Gpara.Pmax./Gpara.adj.*(1-exp(-Gpara.adj.*Gpara.ki.*t_mid./(1+Gpara.adj.*Gpara.kp))./(1+Gpara.adj.*Gpara.kp)),2);
    ki_sys = -(1+adj_sys.*kp_sys)./(adj_sys.*t_mid).*log((1+adj_sys.*sum(u_G_dh.*repmat(Gpara.Pmax./Gpara.adj,dhdata.T,1), 2)).*(1-adj_sys.*innerpara2));   %等效积分系数
    kG_sys = 2*sum(u_G_dh.*Gpara.Pmax.*Gpara.H,2);      %等效惯量
    
    % 电化学储能等效参数求解
    kv_sys = ESpara.Pmax*ESpara.kv;  %等效虚拟惯量系数
    kd_sys = ESpara.Pmax*ESpara.kd;  %等效频率下垂系数
    % 参数求解
    I_sys = kG_sys+kv_sys;  %系统总惯量
    p_sys = adj_sys.*ki_sys./(1+adj_sys.*kp_sys);
    a_sys = 0.5*(I_sys.*adj_sys.*ki_sys+(k_load*dhdata.load+kd_sys).*(1+adj_sys.*kp_sys)+kp_sys)./I_sys./(1+adj_sys.*kp_sys);
    b_sys = ((k_load*dhdata.load+kd_sys).*adj_sys.*ki_sys+ki_sys)./I_sys./(1+adj_sys.*kp_sys);
    
    % 风电等效参数
    kWv_sys = 20*1145.39; % *wp.pmax
    kWd_sys = 12*1145.39; % *wp.pmax
    
    % 系统频差表达式参数求解
    % 系统超调频差(正向功率扰动)
    % t_1 = 14
    % t<=t_1
    rho = kG_sys+kv_sys+kWv_sys; 
    lambda = k_load*dhdata.load+kd_sys+kWd_sys+kp_sys./(1+adj_sys.*kp_sys);
    alpha = adj_sys.*(kG_sys+kv_sys+kWv_sys).^2;
    beta = (1+adj_sys.*kp_sys)./(adj_sys.*ki_sys);
    delta_f = zeros(1, 96);
    omega_d = sqrt((lambda.*beta-rho).^2-4*beta.*rho./alpha)./(2*rho.*beta);
    sigma_1 = (rho+lambda.*beta)./(2*rho.*beta);
    
    rho_1 = kG_sys+kv_sys;
    lambda_1 = k_load*dhdata.load+kd_sys+kp_sys./(1+adj_sys.*kp_sys);
    alpha_1 = adj_sys.*(kG_sys+kv_sys).^2;
    beta_1 = (1+adj_sys.*kp_sys)./(adj_sys.*ki_sys);
    omega_d1 = sqrt((lambda_1.*beta_1-rho_1).^2+4*beta_1.*rho_1./alpha_1)./(2*rho_1.*beta_1);
    sigma_11 = (rho_1+lambda_1.*beta_1)./(2*rho_1.*beta_1);

    for id = 1:dhdata.T
        R_1=(dhdata.el(id)+sum(dhdata.ew(id,:)))*1./(lambda(id)+1./alpha(id));
        R_2=0.5*(dhdata.el(id)+sum(dhdata.ew(id,:)))*(rho(id)-lambda(id).*beta(id)-2*beta(id)./alpha(id))/((lambda(id)+1/alpha(id))*omega_d(id)*rho(id)*beta(id));
        f_t1 = 0;

        for t = 1:0.05:t_1
            t_nadir = 1./omega_d(id)*atanh((2*rho(id)*beta(id)*omega_d(id))./(lambda(id)*beta(id)-rho(id)));
            f_t1 = (R_1.*(cosh(omega_d(id)*t)*exp(-sigma_1(id)*t)-1)+R_2*sinh(omega_d(id)*t)*exp(-sigma_1(id)*t));
            if(abs(f_t1)>abs(delta_f(id)))
                delta_f(id)=f_t1;
            end
            display("T: "+num2str(id)+" t: "+num2str(t)+" df:"+num2str(f_t1));
        end

        for t = t_1:0.05:15
            delta_g1=sigma_1(id)./omega_d(id)*(dhdata.el(id)+sum(dhdata.ew(id,:)))./(lambda(id)+1./(alpha(id)))*sinh(omega_d(id)*t_1)*exp(-sigma_1(id)*t_1)...
                +(dhdata.el(id)+sum(dhdata.ew(id,:)))./(lambda(id)+1./(alpha(id)))*(cosh(omega_d(id)*t_1)*exp(-sigma_1(id)*t_1)-1);
            S_1=f_t1+(dhdata.el(id)+sum(dhdata.ew(id,:))+my_Phi)./(lambda_1(id)+1./(alpha_1(id)));
            S_2=(rho_1(id)-lambda_1(id)*beta_1(id))*f_t1./(2*rho_1(id)*beta_1(id)*omega_d(id))-delta_g1./(alpha_1(id)*rho_1(id)*omega_d1(id))+...
                0.5*(dhdata.el(id)+sum(dhdata.ew(id,:))+my_Phi)*(rho_1(id)-lambda_1(id)*beta_1(id))./((lambda_1(id)+1./(alpha_1(id)))*omega_d1(id)*rho_1(id)*beta_1(id));
            f_t2=exp(sigma_11(id)*(t_1-t))*(S_1*cosh(omega_d1(id)*(t-t_1))+S_2*sinh(omega_d1(id)*(t-t_1)))-(dhdata.el(id)+sum(dhdata.ew(id,:))+my_Phi)./(lambda_1(id)+1./alpha_1(id));
            if(abs(f_t2)>abs(delta_f(id)))
                delta_f(id)=f_t2;
            end
            display("T: "+num2str(id)+" t: "+num2str(t)+" df:"+num2str(f_t2));
        end

    end

    display(delta_f);
    % 系统超调抑制能力下限更新
    for t = 1:dhdata.T
        if abs(delta_f(t))>df.max_delta_f
            display("overshoot exceed standard at t:" + num2str(t));
            Y_syslim(t) = sum(u_G_dh(t,:).*Gpara.Y)+1;
        end
    end
    % 收敛判断
    if sum(abs(delta_f)>df.max_delta_f)==0
        break
    end
end
toc
%% 结果计算及保存(超调频差已在迭代过程完成计算)
df_steady = -0.5*(dhdata.el+sum(dhdata.ew,2)+my_Phi)./(lambda_1+1./alpha_1);    %系统稳态频差(正向功率扰动)
vf_max = -(dhdata.el+sum(dhdata.ew))./rho;                                  %系统最大频率变化率(正向功率扰动)
P_G_PFM = -(1+Gpara.adj.*Gpara.Pmax)./Gpara.kp.*df_steady;                  %火电机组一次调频稳态功率变化量
P_ES_PFM = -ESpara.kv*df_steady;
delta_f=delta_f';
%电化学储能一次调频稳态功率变化量
writematrix(P_G_dh,'dayahead_result.xlsx','sheet','G_result','range','B3:I98')
writematrix(u_G_dh,'dayahead_result.xlsx','sheet','G_result','range','J3:Q98')
writematrix(P_G_PFM,'dayahead_result.xlsx','sheet','G_result','range','R3:Y98')
writematrix(P_ES_dh,'dayahead_result.xlsx','sheet','ES_result','range','B2:B98')
writematrix(E_ES_dh,'dayahead_result.xlsx','sheet','ES_result','range','C2:C98')
writematrix(u_ES_dh,'dayahead_result.xlsx','sheet','ES_result','range','D2:D98')
writematrix(P_ES_PFM,'dayahead_result.xlsx','sheet','ES_result','range','E2:E98')
writematrix(P_wind_dh,'dayahead_result.xlsx','sheet','W_result','range','D3:E98')
writematrix(df_steady,'dayahead_result.xlsx','sheet','f_result','range','B2:B98')
writematrix(delta_f,'dayahead_result.xlsx','sheet','f_result','range','D2:D98')
writematrix(vf_max,'dayahead_result.xlsx','sheet','f_result','range','F2:F98')
writematrix(C_G_run_dh,'dayahead_result.xlsx','sheet','Cost_result','range','A2:A2')
writematrix(C_G_on_dh,'dayahead_result.xlsx','sheet','Cost_result','range','B2:B2')
writematrix(C_wind_dh,'dayahead_result.xlsx','sheet','Cost_result','range','C2:C2')