tic

Date='20190401';
Iterates=10000;
M=200;
a_1=1;
a_2=0.5;
Data=table2array(readtable(strcat(Date,'.csv'),"Format","auto"));
% Data=table2array(readtable("TestData.csv","Format","auto"));

K=size(Data,2)-1;
n=sum(Data(:,end));

Z_Tb=zeros(Iterates,sum(Data(:,end)));
Theta_Tb=cell(size(Data,2)-1,M);
cnt=1;
while cnt<=size(Data,2)-1
    cnt1=1;
    while cnt1<=M
        Carry=zeros(Iterates,1);
        Carry(1,1)=1/2;
        Theta_Tb{cnt,cnt1}=Carry;
        cnt1=cnt1+1;
    end
    cnt=cnt+1;
end
Pi_Tb=zeros(Iterates,M);
Alpha_Tb=zeros(Iterates,1);
N_Tb=zeros(Iterates,1);
Omega_Tb=zeros(Iterates,M);

%LOAD START STATES
Alpha_Tb(1,1)=a_1*a_2;
cnt=2;
Pi_Tb(1,1)=1/(1+a_1*a_2);
while cnt<=M-1
    Pi_Tb(1,cnt)=(1/(1+a_1*a_2))*(1-1/(1+a_1*a_2))^(cnt-1);
    cnt=cnt+1;
end
Pi_Tb(1,M)=(1-1/(1+a_1*a_2))^(M-1);
Z_Tb(1,:)=randsample(M,n,'true',Pi_Tb(1,:));%%%
N_Tb(1,1)=80706;
Omega_Tb(1,:)=round(Pi_Tb(1,:).*(4*n));



cnt=2;
while cnt<=Iterates
    Progress=fprintf('Progress: %.2f%%',cnt*100/Iterates);
    %SAMPLE FROM Z
    cnt1=1;
    cnt4=1;
    while cnt1<=size(Data,1)
        x=Data(cnt1,1:K);
        Prob=zeros(1,M);
        cnt2=1;
        while cnt2<=M
            Prob(cnt2)=Pi_Tb(cnt-1,cnt2);
            cnt3=1;
            while cnt3<=K
                Prob(cnt2)=Prob(cnt2)*Theta_Tb{cnt3,cnt2}(cnt-1)^...
                        (x(cnt3))*(1-Theta_Tb{cnt3,cnt2}(cnt-1))^(1-x(cnt3));
                cnt3=cnt3+1;
            end
            cnt2=cnt2+1;
        end
        Prob=Prob./sum(Prob);
        Z_Tb(cnt,cnt4:cnt4+Data(cnt1,end)-1)=randsample(M,Data(cnt1,end),'true',Prob)';
        cnt4=cnt4+Data(cnt1,end);
        cnt1=cnt1+1;
    end

    %SAMPLE FROM THETA
    n_k=zeros(1,M);
    n_jk=zeros(K,M);
    cnt1=1;
    while cnt1<=n
        n_k(Z_Tb(cnt,cnt1))=n_k(Z_Tb(cnt,cnt1))+1;
        x=Data(useful(cnt1,Data),1:size(Data,2)-1);
        cnt2=1;
        while cnt2<=K
            if x(cnt2)==1
                n_jk(cnt2,Z_Tb(cnt,cnt1))=n_jk(cnt2,Z_Tb(cnt,cnt1))+1;
            end
            cnt2=cnt2+1;
        end
        cnt1=cnt1+1;
    end
    cnt1=1;
    while cnt1<=K
        cnt2=1;
        while cnt2<=M
            Theta_Tb{cnt1,cnt2}(cnt)=betarnd(n_jk(cnt1,cnt2)+1,...
                n_k(cnt2)-n_jk(cnt1,cnt2)+Omega_Tb(cnt-1,cnt2)+1);
            cnt2=cnt2+1;
        end
        cnt1=cnt1+1;
    end

    %SAMPLE FROM PI
    c_k=n_k+Omega_Tb(cnt-1,:);
    cnt1=1;
    Carry=zeros(1,M);
    while cnt1<=M-1
        Carry(cnt1)=betarnd(1+c_k(cnt1),Alpha_Tb(cnt-1,1)+sum(c_k(cnt1+1:M)));
        cnt1=cnt1+1;
    end
    Carry(end)=1;
    cnt1=1;
    Pi_Tb(cnt,:)=ones(1,M);
    while cnt1<=M-1
        Pi_Tb(cnt,cnt1)=Pi_Tb(cnt,cnt1)*Carry(cnt1);
        Pi_Tb(cnt,cnt1+1:end)=Pi_Tb(cnt,cnt1+1:end).*(1-Carry(cnt1));
        cnt1=cnt1+1;
    end

    %SAMPLE FROM ALPHA
    Alpha_Tb(cnt,1)=gamrnd(a_1-1+M,a_2-log(Pi_Tb(cnt,M)));

    %SAMPLE FROM N
    Carry=zeros(1,M);
    cnt1=1;
    while cnt1<=M
        Carry(cnt1)=Pi_Tb(cnt,cnt1);
        cnt2=1;
        while cnt2<=K
            Carry(cnt1)=Carry(cnt1)*(1-Theta_Tb{cnt2,cnt1}(cnt));
            cnt2=cnt2+1;
        end
        cnt1=cnt1+1;
    end
    N_Tb(cnt)=nbinrnd(n,1-sum(Carry))+n;

    %SAMPLE FROM OMEGA
    Omega_Tb(cnt,:)=round(mnrnd(N_Tb(cnt)-n,Carry/sum(Carry)));
    cnt=cnt+1;
    fprintf(repmat('\b',1,Progress))
end

Data=N_Tb(Iterates-round(Iterates/2):Iterates);
pd=fitdist(Data,'kernel','Kernel','Normal');
disp('_______________________')
disp(strcat('Data used: ',Date))
fprintf('Iterates = %u \n',Iterates)
fprintf('M = %u \n',M)

M=max(Data);
m=min(Data);
X=m:(M-m)/999:M;
x=5*m/4-M/4:(3/2)*(M-m)/9999:5*M/4-m/4;
y=pdf(pd,x);
Y=pdf(pd,X);
figure
plot(x,y,'Color','Black')
xlabel("N")
figure
plot(X,Y,'Color','Black')
xlabel("N")

fprintf('Min = %u \n',m)

Quant=cumtrapz(x,y);
QuantCheck=[0.05,0.5,0.95];
cnt=1;
while cnt<=size(Quant,2)-1
    cnt1=1;
    while cnt1<=size(QuantCheck,2)
        if and(Quant(cnt)<=QuantCheck(cnt1),Quant(cnt+1)>=QuantCheck(cnt1))==1
            fprintf('%.2f Quantile = %.1f \n',QuantCheck(cnt1),...
                round((cnt+1/2)/size(Quant,2)*(3/2*(M-m))+5*m/4-M/4,1))
        end
        cnt1=cnt1+1;
    end
    cnt=cnt+1;
end
fprintf('Max = %u \n',M)
fprintf('Mean = %.1f \n',mean(pd))
fprintf('Std = %.1f \n',std(pd))
fprintf('Observed cases n = %u \n',n)
toc
disp('_______________________')
beep on
beep

function counter = useful(number,Data)
sum=0;
counter=0;
while sum<number
    sum=sum+Data(counter+1,end);
    counter=counter+1;
end
end