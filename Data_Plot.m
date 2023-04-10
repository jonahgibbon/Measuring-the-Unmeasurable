Smth=1;
XTime=[
    datetime("2018-04-01")
    datetime("2018-07-01")
    datetime("2018-10-01")
    datetime("2019-01-01")
    datetime("2019-04-01")
    datetime("2019-07-01")
    datetime("2019-10-01")    
    datetime("2020-01-01")
    datetime("2020-04-01")
    datetime("2020-07-01")
    datetime("2020-10-01")
    datetime("2021-01-01")
    ];
Data=[
    32312 34552.0 36163.4 37258.8 38326.4 40123.3 42612   
    27754 30163.5 31559.6 32517.5 33529.6 35055.7 37271
    36079 37960.6 39563.4 40770.2 42073.2 44090.8 48650
    34939 39115.9 41118.8 42514.8 44013.7 46228.9 49235
    32761 36348.3 38140.2 39465.5 40683.2 42509.8 44914
    34064 36700.9 38473.7 39818.2 41126.3 43175.7 46869
    25665 26841.4 27830.8 28476.7 29245.5 30650.4 33051
    19167 20125.2 20779.7 21277.7 21805.0 22733.3 24182
    26411 27972.4 29003.7 29812.1 30667.4 31982.1 33356
    34222 36178.4 37540.9 38577.6 39615.8 41255.6 44326
    23731 24403.3 25074.0 25592.0 26159.9 27130.6 29093
    24599 25814.6 26570.5 27105.2 27689.8 28615.4 30339
    ];
XNum=zeros(size(XTime,1),1);
cnt=1;
while cnt<=size(XTime,1)
    XNum(cnt)=days(XTime(cnt,1)-XTime(1,1));
    cnt=cnt+1;
end
Cve_Mi=csaps(XNum,Data(:,1)',Smth);
Cve_05=csaps(XNum,Data(:,2)',Smth);
Cve_25=csaps(XNum,Data(:,3)',Smth);
Cve_50=csaps(XNum,Data(:,4)',Smth);
Cve_75=csaps(XNum,Data(:,5)',Smth);
Cve_95=csaps(XNum,Data(:,6)',Smth);
Cve_Ma=csaps(XNum,Data(:,7)',Smth);
Cve_Stright=csaps(XNum,Data(:,4)',0);


figure
scatter(XTime,Data(:,4),'ko','filled','LineWidth',1)
hold on
plot(linspace(0,1006,1000),ppval(Cve_05,linspace(0,1006,1000)),'k--','LineWidth',0.75)
plot(linspace(0,1006,1000),ppval(Cve_50,linspace(0,1006,1000)),'k-','LineWidth',1)
plot(linspace(0,1006,1000),ppval(Cve_95,linspace(0,1006,1000)),'k--','LineWidth',0.75)
plot(linspace(0,1006,1000),ppval(Cve_Stright,linspace(0,1006,1000)),'-','LineWidth',0.5,'Color',[.7 .7 .7])
xlabel("Date")
ylabel("N")
ylim([1.8*10^4,5*10^4])
legend(["","5% Quantile","50% Quantile","95% Quantile"],'Location','northeast')
print('Results','-depsc')

latex(sym(vpa(Data)))
