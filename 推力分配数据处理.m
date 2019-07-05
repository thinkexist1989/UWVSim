%% 推力分配数据处理

%% 姿态数据加噪声
position = Eta_Eth.Data(:,1:3)*1000+10*randn(251,3);
attitude = Eta_Eth.Data(:,4:6)*180/pi+1*randn(251,3);
