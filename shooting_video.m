%%%%%%%%%%%%%%%%%%%%%
% シューティング法
% 日野将志
%%%%%%%%%%%%%%%%%%%%%

clear;clc
%% パラメータ
alpha=0.33;% 資本シェア
beta=0.96;% 時間割引
delta=0.08;% 資本減耗
T=50;% シミュレーション期間
tol=1e-4;%収束の許容できる誤差

%% 定常状態値
kss=alpha/(1/beta-(1-delta));
kss=kss^(1/(1-alpha));%資本の定常状態値
css=kss^alpha-delta*kss;%消費の定常状態値

%% 初期値

k0=0.2*kss;%資本の初期値
c0l=tol;%消費の初期値の下限
c0u=css;%消費の初期値の上限

%% solving model with shooting
for i=1:100
    kseq=zeros(T,1);
    cseq=kseq;
    kseq(1)=k0;
    
    cseq0=0.5*(c0l+c0u);
    cseq(1)=cseq0;
    for t=1:T-1
        if cseq(t)<c0l-tol
            %c<0にならないように．位相図で右下に発散に相当
            cseq(t)=0;
            break
        end
        if kseq(t)<k0
            %k<0にならないように．位相図で左上に発散に相当
            kseq(t)=0;
            break
        end
        kseq(t+1)=kseq(t)^alpha+(1-delta)*kseq(t)-cseq(t);
        cseq(t+1)=beta*cseq(t)*(1+alpha*kseq(t+1)^(alpha-1)-delta);
    end
    if cseq(t)-css>tol
        c0u=cseq0;
    else
        c0l=cseq0;
    end
    
    if abs(kseq(T)-kss)<tol
        break
    end
    display1=sprintf('%4d th iteration',i);
    disp(display1)
    fprintf('kseq(T), kss, cseq(T), cseq0,\n')
    display=[kseq(T),kss,cseq(T),cseq0];
    disp(display)
end

%% divergent path

cseq1=zeros(T,1);
kseq1=zeros(T,1);
cseq2=zeros(T,1);
kseq2=zeros(T,1);

cseq1(1)=cseq(1)-0.005;
kseq1(1)=kseq(1);
for t=1:T-1
    if cseq1(t)<0
        %c<0にならないように．位相図で右下に発散に相当
        cseq(t)=0;
        break
    end
    if kseq1(t)<0
        %k<0にならないように．位相図で左上に発散に相当
        kseq(t)=0;
        break
    end
    kseq1(t+1)=kseq1(t)^alpha+(1-delta)*kseq1(t)-cseq1(t);
    cseq1(t+1)=beta*cseq1(t)*(1+alpha*kseq1(t+1)^(alpha-1)-delta);
end

cseq2(1)=cseq(1)+0.00004;
kseq2(1)=kseq(1);
for t=1:T-1
    if cseq1(t)<0
        %c<0にならないように．位相図で右下に発散に相当
        cseq(t)=0;
        break
    end
    if kseq1(t)<0
        %k<0にならないように．位相図で左上に発散に相当
        kseq(t)=0;
        break
    end
    kseq2(t+1)=kseq2(t)^alpha+(1-delta)*kseq2(t)-cseq2(t);
    cseq2(t+1)=beta*cseq2(t)*(1+alpha*kseq2(t+1)^(alpha-1)-delta);
end


%% static plot

kdomain=[0:0.01:9.2*kss];
index=sum(kdomain<kss);
Deltac0=kdomain(index);
Deltak0=kdomain.^alpha-delta.*kdomain;

figure
plot([kss,kss],[0,3],'-.',kdomain,Deltak0,'-.')
hold on
plot(kseq,cseq,'.',kseq1,cseq1,'.',kseq2,cseq2,'.')
legend('\Delta c=0','\Delta k=0','Dynamics','Divergence1','Divergence2')

%% video part
video = VideoWriter('Convergence');
open(video);

figure
plot([kss,kss],[0,3],'-.',kdomain,Deltak0,'-.')
hold on
for t=1:T-1
    kmovie= kseq(t);
    cmovie=cseq(t);
    plot(kmovie,cmovie,'oblack','LineWidth',0.3)
    hold on
    %u.Value = i;
    frame=getframe(gcf);
    writeVideo(video,frame)
end
legend('\Delta c=0','\Delta k=0','Dynamics','Divergence1','Divergence2')
title('convergence','Fontsize',12)
xlabel('capital','Fontsize',12)
%hold off
close(video)


