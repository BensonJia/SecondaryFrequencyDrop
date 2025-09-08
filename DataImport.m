%数据导入
%% 火电机组参数
Gpara.N = 8;                                                                %火电机组数量
Gpara.Pmax = readmatrix('parameter.xlsx','Sheet','unit','Range','E3:L3');   %额定出力
Gpara.Pmin = readmatrix('parameter.xlsx','Sheet','unit','Range','E4:L4');   %最小技术出力
Gpara.on = readmatrix('parameter.xlsx','Sheet','unit','Range','E5:L5');     %开机成本
Gpara.sst = readmatrix('parameter.xlsx','Sheet','unit','Range','E6:L6');    %最小启停机时间
Gpara.ramp = readmatrix('parameter.xlsx','Sheet','unit','Range','E7:L7');   %爬坡率
Gpara.a = readmatrix('parameter.xlsx','Sheet','unit','Range','E8:L8');      %购电成本系数a
Gpara.b = readmatrix('parameter.xlsx','Sheet','unit','Range','E9:L9');      %购电成本系数b
Gpara.c = readmatrix('parameter.xlsx','Sheet','unit','Range','E10:L10');    %负荷正备用成本系数
Gpara.d = readmatrix('parameter.xlsx','Sheet','unit','Range','E11:L11');    %负荷负备用成本系数
Gpara.adj = readmatrix('parameter.xlsx','Sheet','unit','Range','E12:L12');  %调差系数
Gpara.kp = readmatrix('parameter.xlsx','Sheet','unit','Range','E13:L13');   %调速器比例增益
Gpara.ki = readmatrix('parameter.xlsx','Sheet','unit','Range','E14:L14');   %调速器积分增益
Gpara.H = readmatrix('parameter.xlsx','Sheet','unit','Range','E15:L15');    %惯性时间常数
Gpara.Y = readmatrix('parameter.xlsx','Sheet','unit','Range','E16:L16');    %超调抑制能力
Gpara.delay = readmatrix('parameter.xlsx','Sheet','unit','Range','E17:L17');%事故备用启动时延
Gpara.stab = readmatrix('parameter.xlsx','Sheet','unit','Range','E18:L18'); %调速器稳定时间
%% 电化学储能参数
ESpara.N = 1;                                                               %电化学储能数量
ESpara.Pmax = readmatrix('parameter.xlsx','Sheet','unit','Range','O7:O7');  %额定功率
ESpara.Emax = readmatrix('parameter.xlsx','Sheet','unit','Range','O8:O8');  %额定容量
ESpara.Emin = readmatrix('parameter.xlsx','Sheet','unit','Range','O9:O9');  %能量下限
ESpara.Eini = readmatrix('parameter.xlsx','Sheet','unit','Range','O10:O10');%初始能量
ESpara.eff = readmatrix('parameter.xlsx','Sheet','unit','Range','O11:O11'); %功率转换效率
ESpara.kv = readmatrix('parameter.xlsx','Sheet','unit','Range','O12:O12');  %虚拟惯量系数
ESpara.kd = readmatrix('parameter.xlsx','Sheet','unit','Range','O13:O13');  %频率下垂系数
%% 日前负荷及风电
dhdata.T = 96;                                                                  %日前时段数
dhdata.load = readmatrix('load.xlsx','Sheet','intraday','Range','E2:E97');      %日前负荷预测值
dhdata.wind(:,1) = readmatrix('wind1.xlsx','Sheet','intraday','Range','D2:D97');%日前风电1预测值
dhdata.wind(:,2) = readmatrix('wind2.xlsx','Sheet','intraday','Range','D2:D97');%日前风电2预测值
dhdata.el = readmatrix('load.xlsx','Sheet','intraday','Range','G2:G97');        %日前负荷预测误差
dhdata.ew(:,1) = readmatrix('wind1.xlsx','Sheet','intraday','Range','F2:F97');  %日前风电1预测误差
dhdata.ew(:,2) = readmatrix('wind2.xlsx','Sheet','intraday','Range','F2:F97');  %日前风电2预测误差
%% 日内负荷及风电
for t = 1:24
    iddata(t).T = 96-4*(t-1);                                                   %日内时段数
    iddata(t).load = readmatrix('load.xlsx','Sheet',t+2,'Range','C:C');         %日内负荷预测值
    iddata(t).wind(:,1) = readmatrix('wind1.xlsx','Sheet',t+2,'Range','C:C');   %日内风电1预测值
    iddata(t).wind(:,2) = readmatrix('wind2.xlsx','Sheet',t+2,'Range','C:C');   %日内风电2预测值
    iddata(t).el = readmatrix('load.xlsx','Sheet',t+2,'Range','E:E');           %日内负荷预测误差
    iddata(t).ew(:,1) = readmatrix('wind1.xlsx','Sheet',t+2,'Range','E:E');     %日内风电1预测误差
    iddata(t).ew(:,2) = readmatrix('wind2.xlsx','Sheet',t+2,'Range','E:E');     %日内风电2预测误差
    iddata(t).load(1) = [];
    iddata(t).wind(1,:) = [];
    iddata(t).el(1) = [];
    iddata(t).ew(1,:) = [];
end
%% 其他参数
M = 100000000;          %大M参数
df.steady1 = 0.004;     %仅有负荷与风电波动时的稳态频差上限标幺值(有名值为0.2Hz)
df.overshoot1 = 0.01;   %仅有负荷与风电波动时的超调频差上限标幺值(有名值为0.5Hz)
df.changerate1 = 0.004; %仅有负荷与风电波动时的频率变化率上限标幺值(有名值为0.2Hz/s)
df.steady2 = 0.01;      %单机故障与负荷、风电波动同时发生时的稳态频差上限标幺值(有名值为0.5Hz)
df.overshoot2 = 0.015;  %单机故障与负荷、风电波动同时发生时的超调频差上限标幺值(有名值为0.75Hz)
df.changerate2 = 0.01;  %单机故障与负荷、风电波动同时发生时的频率变化率上限标幺值(有名值为0.5Hz/s)
df.deadzone = 0.00066;  %调频死区标幺值(有名值为0.033Hz)
df.max_delta_f = 1;   %仅有负荷与风电波动时的超调频差上限(0.5Hz)
k_load = 1;             %负荷频率响应系数
penalty_wind = 200;     %弃风成本
N_wind = 2;             %风电场数量
%% 结果保存
clear t
save dataimport