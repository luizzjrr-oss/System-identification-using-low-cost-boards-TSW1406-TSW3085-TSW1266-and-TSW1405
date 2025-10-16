 clc
clear
close all
format shortG
format compact
% format short

% 
% fc1=0/307.2;
% fc2=150/307.2;

figure
hold on

for k_max_o_k=4:1:4
for k_max_d_k=4:1:4
for k_tip_Amp=2:1:2
for k_t_amp_pot=1.00:-0.25:1.00;

    
vet_nmse_t_amp_pot_k=[];
vet_nmse_t_amp_pot_o=[];

k_t_amp_pot_k=0;

for k_Norm_k=1.0:40.0:1.0
for k_Norm_k2=1.0:40.0:1.0
k_t_amp_pot_k=k_t_amp_pot_k+1;    

   
MAToo=[1 0 1;...
       2 0 1;...
       3 0 1;...
       4 2 1;...
       1 2 1;...
       2 0 2;...
       3 2 1;...
       4 2 2;...
       1 0 2;...
       4 0 2;...
       3 2 2;...
       4 0 1;...
       ];    


    
Sinal3BB1=[];Sinal3BB2=[];Sinal3BB3=[];Sinal_r3BB1=[];Sinal_r3BB2=[];Sinal_r3BB3=[];

if k_tip_Amp==1
% 
load(['Signalx1t1.mat']);
Sinal3BB1=[Sinal3BB1;SinalBB3];
Sinal_r3BB1=[Sinal_r3BB1;SinalBBdmod3];
load(['Signalx2t1.mat']);
Sinal3BB1=[Sinal3BB1;SinalBB3];
Sinal_r3BB1=[Sinal_r3BB1;SinalBBdmod3];
load(['Signalx3t1.mat']);
Sinal3BB1=[Sinal3BB1;SinalBB3];
Sinal_r3BB1=[Sinal_r3BB1;SinalBBdmod3];
load(['Signalx4t1.mat']);
Sinal3BB1=[Sinal3BB1;SinalBB3];
Sinal_r3BB1=[Sinal_r3BB1;SinalBBdmod3];
load(['Signalx5t1.mat']);
Sinal3BB1=[Sinal3BB1;SinalBB3];
Sinal_r3BB1=[Sinal_r3BB1;SinalBBdmod3];
load(['Signalx1t2.mat']);
Sinal3BB1=[Sinal3BB1;SinalBB3];
Sinal_r3BB1=[Sinal_r3BB1;SinalBBdmod3];

load(['Signalx1t2.mat']);
Sinal3BB2=[Sinal3BB2;SinalBB3];
Sinal_r3BB2=[Sinal_r3BB2;SinalBBdmod3];
load(['Signalx2t2.mat']);
Sinal3BB2=[Sinal3BB2;SinalBB3];
Sinal_r3BB2=[Sinal_r3BB2;SinalBBdmod3];
load(['Signalx3t2.mat']);
Sinal3BB2=[Sinal3BB2;SinalBB3];
Sinal_r3BB2=[Sinal_r3BB2;SinalBBdmod3];
load(['Signalx4t2.mat']);
Sinal3BB2=[Sinal3BB2;SinalBB3];
Sinal_r3BB2=[Sinal_r3BB2;SinalBBdmod3];
load(['Signalx5t2.mat']);
Sinal3BB2=[Sinal3BB2;SinalBB3];
Sinal_r3BB2=[Sinal_r3BB2;SinalBBdmod3];
load(['Signalx1t1.mat']);
Sinal3BB2=[Sinal3BB2;SinalBB3];
Sinal_r3BB2=[Sinal_r3BB2;SinalBBdmod3];

load(['Signalx1t1.mat']);
Sinal3BB3=[Sinal3BB3;SinalBB3];
Sinal_r3BB3=[Sinal_r3BB3;SinalBBdmod3];
load(['Signalx2t1.mat']);
Sinal3BB3=[Sinal3BB3;SinalBB3];
Sinal_r3BB3=[Sinal_r3BB3;SinalBBdmod3];
load(['Signalx3t1.mat']);
Sinal3BB3=[Sinal3BB3;SinalBB3];
Sinal_r3BB3=[Sinal_r3BB3;SinalBBdmod3];
load(['Signalx4t1.mat']);
Sinal3BB3=[Sinal3BB3;SinalBB3];
Sinal_r3BB3=[Sinal_r3BB3;SinalBBdmod3];
load(['Signalx5t1.mat']);
Sinal3BB3=[Sinal3BB3;SinalBB3];
Sinal_r3BB3=[Sinal_r3BB3;SinalBBdmod3];
load(['Signalx1t2.mat']);
Sinal3BB3=[Sinal3BB3;SinalBB3];
Sinal_r3BB3=[Sinal_r3BB3;SinalBBdmod3];

%%% Aqui

elseif k_tip_Amp==2

load(['a2Signalx1t1.mat']);
Sinal3BB1=[Sinal3BB1;SinalBB3];
Sinal_r3BB1=[Sinal_r3BB1;SinalBBdmod3];
load(['a2Signalx2t1.mat']);
Sinal3BB1=[Sinal3BB1;SinalBB3];
Sinal_r3BB1=[Sinal_r3BB1;SinalBBdmod3];
load(['a2Signalx3t1.mat']);
Sinal3BB1=[Sinal3BB1;SinalBB3];
Sinal_r3BB1=[Sinal_r3BB1;SinalBBdmod3];
load(['a2Signalx4t1.mat']);
Sinal3BB1=[Sinal3BB1;SinalBB3];
Sinal_r3BB1=[Sinal_r3BB1;SinalBBdmod3];
load(['a2Signalx5t1.mat']);
Sinal3BB1=[Sinal3BB1;SinalBB3];
Sinal_r3BB1=[Sinal_r3BB1;SinalBBdmod3];
load(['a2Signalx1t2.mat']);
Sinal3BB1=[Sinal3BB1;SinalBB3];
Sinal_r3BB1=[Sinal_r3BB1;SinalBBdmod3];

load(['a2Signalx1t2.mat']);
Sinal3BB2=[Sinal3BB2;SinalBB3];
Sinal_r3BB2=[Sinal_r3BB2;SinalBBdmod3];
load(['a2Signalx2t2.mat']);
Sinal3BB2=[Sinal3BB2;SinalBB3];
Sinal_r3BB2=[Sinal_r3BB2;SinalBBdmod3];
load(['a2Signalx3t2.mat']);
Sinal3BB2=[Sinal3BB2;SinalBB3];
Sinal_r3BB2=[Sinal_r3BB2;SinalBBdmod3];
load(['a2Signalx4t2.mat']);
Sinal3BB2=[Sinal3BB2;SinalBB3];
Sinal_r3BB2=[Sinal_r3BB2;SinalBBdmod3];
load(['a2Signalx5t2.mat']);
Sinal3BB2=[Sinal3BB2;SinalBB3];
Sinal_r3BB2=[Sinal_r3BB2;SinalBBdmod3];
load(['a2Signalx1t1.mat']);
Sinal3BB2=[Sinal3BB2;SinalBB3];
Sinal_r3BB2=[Sinal_r3BB2;SinalBBdmod3];

load(['a2Signalx1t1.mat']);
Sinal3BB3=[Sinal3BB3;SinalBB3];
Sinal_r3BB3=[Sinal_r3BB3;SinalBBdmod3];
load(['a2Signalx2t1.mat']);
Sinal3BB3=[Sinal3BB3;SinalBB3];
Sinal_r3BB3=[Sinal_r3BB3;SinalBBdmod3];
load(['a2Signalx3t1.mat']);
Sinal3BB3=[Sinal3BB3;SinalBB3];
Sinal_r3BB3=[Sinal_r3BB3;SinalBBdmod3];
load(['a2Signalx4t1.mat']);
Sinal3BB3=[Sinal3BB3;SinalBB3];
Sinal_r3BB3=[Sinal_r3BB3;SinalBBdmod3];
load(['a2Signalx5t1.mat']);
Sinal3BB3=[Sinal3BB3;SinalBB3];
Sinal_r3BB3=[Sinal_r3BB3;SinalBBdmod3];
load(['a2Signalx1t2.mat']);
Sinal3BB3=[Sinal3BB3;SinalBB3];
Sinal_r3BB3=[Sinal_r3BB3;SinalBBdmod3];

end

passo_delay=1;

ini_t_delay=k_max_d_k;%%%% Delay mínimo do algoritimo
max_delay=k_max_d_k;%%%% Delay máxima do algoritimo
lim_casas=8;%%%% limite de casas dos coeficientes
N_P=length(Sinal3BB1);%%%% Número de pontos dos dados
ini_nmse=1000;%%%% Corta o inicio dos sinais
ini_u_ordem=k_max_o_k;%%%% Ordem mú‹ima do algoritimo
max_u_ordem=k_max_o_k;%%%% Ordem máxima do algoritimo
nomali_n=1;%%%% 1== sinal normalizado 0= não normalizado
retira_dc=0;%%%%1== retira componente DC
bandaS=[0.00 1.00];%%% Banda dos sinais (Classificação, Validação e Sinal) Não utilizar limites de bandas menores que 0.15 e maiores que 0.85
ordem_C=7;%%% Ordem do corte dos sinais acima
bandaS2=[0.00 1.00];%%% Banda dos sinais (Classificação, Validação e Sinal) Não utilizar limites de bandas menores que 0.15 e maiores que 0.85
ordem_C2=25;%%% Ordem do corte dos sinais acima
k_mul_a=1;%%% Constante pela qual todos os valores serão multiplicados
k_amp_s=1;%%% Amplitude da simulação
Acha_corr=0;%%% 1== sincroniza os dois sinais na máxima correlação 0= não sincroniza
sim_atr=0;%%%% Força um atraso no sinal (simulação)
procura_atr_fino=0;%%%% Procura um atraso fino
sem_DC=1;
sem_REALIM=1;
passo_int=1;
k_IQ=1;%%% = (1,IQ) (=2,I) (=3,Q)
iKarnelV=1;
fKarnelV=1;

S_DD_Fit=real(Sinal_r3BB1)';
S_DD_Fiv=real(Sinal_r3BB2)';
S_DD_Fis=real(Sinal_r3BB3)';

Q_DD_Fit=imag(Sinal_r3BB1)';
Q_DD_Fiv=imag(Sinal_r3BB2)';
Q_DD_Fis=imag(Sinal_r3BB3)';

S_A_Fit=real(Sinal3BB1)';
S_A_Fiv=real(Sinal3BB2)';
S_A_Fis=real(Sinal3BB3)';

Q_A_Fit=imag(Sinal3BB1)';
Q_A_Fiv=imag(Sinal3BB2)';
Q_A_Fis=imag(Sinal3BB3)';

% S_DD_Fit=real(Sinal3BB1)';
% S_DD_Fiv=real(Sinal3BB2)';
% S_DD_Fis=real(Sinal3BB3)';
% 
% Q_DD_Fit=imag(Sinal3BB1)';
% Q_DD_Fiv=imag(Sinal3BB2)';
% Q_DD_Fis=imag(Sinal3BB3)';
% 
% S_A_Fit=real(Sinal_r3BB1)';
% S_A_Fiv=real(Sinal_r3BB2)';
% S_A_Fis=real(Sinal_r3BB3)';
% 
% Q_A_Fit=imag(Sinal_r3BB1)';
% Q_A_Fiv=imag(Sinal_r3BB2)';
% Q_A_Fis=imag(Sinal_r3BB3)';

% 
% 

a_fir1 = [];
a_fir2 = [];
a_fir3 = [];
a_fir4 = [];
a_fir5 = [];
a_fir6 = [];
a_fir7 = [];
a_fir8 = [];
b_fir = [];
lib_fir = [];
Sx_fir=[];
Sy_fir=[];

if k_IQ==1

S_DD_Fit=transpose(S_DD_Fit)+transpose(Q_DD_Fit*1i);
S_DD_Fiv=transpose(S_DD_Fiv)+transpose(Q_DD_Fiv*1i);
S_DD_Fis=transpose(S_DD_Fis)+transpose(Q_DD_Fis*1i);

S_A_Fit=transpose(S_A_Fit)+transpose(Q_A_Fit*1i);
S_A_Fiv=transpose(S_A_Fiv)+transpose(Q_A_Fiv*1i);
S_A_Fis=transpose(S_A_Fis)+transpose(Q_A_Fis*1i);

end
if k_IQ==2

S_DD_Fit=transpose(S_DD_Fit);%+transpose(Q_DD_Fit*1i);
S_DD_Fiv=transpose(S_DD_Fiv);%+transpose(Q_DD_Fiv*1i);
S_DD_Fis=transpose(S_DD_Fis);%+transpose(Q_DD_Fis*1i);

S_A_Fit=transpose(S_A_Fit);%+transpose(Q_A_Fit*1i);
S_A_Fiv=transpose(S_A_Fiv);%+transpose(Q_A_Fiv*1i);
S_A_Fis=transpose(S_A_Fis);%+transpose(Q_A_Fis*1i);

end
if k_IQ==3
    
S_DD_Fit=transpose(Q_DD_Fit);%+transpose(Q_DD_Fit*1i);
S_DD_Fiv=transpose(Q_DD_Fiv);%+transpose(Q_DD_Fiv*1i);
S_DD_Fis=transpose(Q_DD_Fis);%+transpose(Q_DD_Fis*1i);

S_A_Fit=transpose(Q_A_Fit);%+transpose(Q_A_Fit*1i);
S_A_Fiv=transpose(Q_A_Fiv);%+transpose(Q_A_Fiv*1i);
S_A_Fis=transpose(Q_A_Fis);%+transpose(Q_A_Fis*1i);

end
nn = 800;

A=1;

if bandaS(1)==0 && bandaS(2)==1
B=1;
elseif bandaS(1)==0 && bandaS(2)<1
% [B,A] = BUTTER(ordem_C,bandaS(2),'low');    
ff = [0 bandaS(2)-0.02 bandaS(2)-0.01 bandaS(2) bandaS(2)+0.01 1];
aa = [1 1 1 0 0];
up = [1.001 1.001 1.001 0.001 0.001];
lo = [0.999 0.999 0.999 -0.001 -0.001];
B = fircls(nn,ff,aa,up,lo);

elseif bandaS(1)>0 && bandaS(2)==1
% [B,A] = BUTTER(ordem_C,bandaS(1),'high');
ff = [0 bandaS(1)-0.01 bandaS(1) bandaS(1)+0.01 bandaS(1)+0.02 1];
aa = [0 0 1 1 1];
up = [0.001 0.001 1.001 1.001 1.001];
lo = [-0.001 -0.001 0.999 0.999 0.999];
B = fircls(nn,ff,aa,up,lo);

elseif bandaS(1)>0 && bandaS(2)<1
% [B,A] = BUTTER(ordem_C,bandaS);  
ff = [0 bandaS(1)-0.01 bandaS(1) bandaS(2) bandaS(2)+0.01 1];
aa = [0 0 1 0 0];
up = [0.001 0.001 1.001 0.001 0.001];
lo = [-0.001 -0.001 0.999 -0.001 -0.001];
B = fircls(nn,ff,aa,up,lo);

end

S_DD_Fit = filter(B,A,S_DD_Fit);
S_DD_Fiv = filter(B,A,S_DD_Fiv);
S_DD_Fis = filter(B,A,S_DD_Fis);
S_A_Fit = filter(B,A,S_A_Fit);
S_A_Fiv = filter(B,A,S_A_Fiv);
S_A_Fis = filter(B,A,S_A_Fis);

% figure
% hold on
% plot(0:2/length(S_A_Fit):2-2/length(S_A_Fit),db(abs(fft(S_A_Fit))),'b')
% plot(0:2/length(S_DD_Fit):2-2/length(S_DD_Fit),db(abs(fft(S_DD_Fit))),'r')


S_A_Fit = k_mul_a*S_A_Fit;
S_DD_Fit= k_mul_a*S_DD_Fit;
S_A_Fiv = k_mul_a*S_A_Fiv;
S_DD_Fiv= k_mul_a*S_DD_Fiv;
S_A_Fis = k_mul_a*S_A_Fis;
S_DD_Fis= k_mul_a*S_DD_Fis;

if retira_dc==1
S_DD_Fit=S_DD_Fit-mean(S_DD_Fit);
S_A_Fit=S_A_Fit-mean(S_A_Fit);

S_DD_Fiv=S_DD_Fiv-mean(S_DD_Fiv);
S_A_Fiv=S_A_Fiv-mean(S_A_Fiv);

S_DD_Fis=S_DD_Fis-mean(S_DD_Fis);
S_A_Fis=S_A_Fis-mean(S_A_Fis);
end

if nomali_n==1

NormA=sum(abs(S_A_Fit))/(10^4);
NormDD=sum(abs(S_DD_Fit))/(10^4);

NormA=NormA*k_Norm_k;
NormDD=NormDD*k_Norm_k2;

S_A_Fit=S_A_Fit/NormA;
S_DD_Fit=S_DD_Fit/NormDD;

S_A_Fiv=S_A_Fiv/NormA;
S_DD_Fiv=S_DD_Fiv/NormDD;

S_A_Fis=S_A_Fis/NormA;
S_DD_Fis=S_DD_Fis/NormDD;

end
% b3=[1 0];
% a3=[1 0 0]/(b3(1));
% Sx=S_A_Fit;
% for k=max([length(a3) length(b3)]):1:length(Sx)
%    Sy(k)=0;
%     for ka=1:1:length(a3)
%         Sy(k)=Sy(k)+Sx(k-ka+1)*a3(ka);
%     end
%     for kb=2:1:length(b3) %começa do segundo termo 
%         Sy(k)=Sy(k)-Sy(k-kb+1)*b3(kb);
%     end
% end
% S_D_Ce=Sy';
% Desloc=-1;

if Acha_corr==1

[y_cor,x_cor]=max(xcorr(S_A_Fit,S_DD_Fit));
Tam_S=length(S_A_Fit);
Desloc=length(S_A_Fit)-x_cor;
% Desloc=8;
Final_S=Tam_S-Desloc-1;
inici_S=1+Desloc;
S_DD_Fit=S_DD_Fit(inici_S+abs(Desloc):-2-abs(Desloc)+N_P);
S_A_Fit=S_A_Fit(1+abs(Desloc):-abs(Desloc)+Final_S-1);
S_DD_Fiv=S_DD_Fiv(inici_S+abs(Desloc):-2-abs(Desloc)+N_P);
S_A_Fiv=S_A_Fiv(1+abs(Desloc):-abs(Desloc)+Final_S-1);
S_DD_Fis=S_DD_Fis(inici_S+abs(Desloc):-2-abs(Desloc)+N_P);
S_A_Fis=S_A_Fis(1+abs(Desloc):-abs(Desloc)+Final_S-1);
end

% figure 
% hold on
% plot(real(S_DD_Fit),'r')
% plot(real(S_A_Fit),'k')
% legend('entrada A/D','Saida FPGA')
% title('treinamento real')
% figure 
% hold on
% plot(imag(S_DD_Fit),'r')
% plot(imag(S_A_Fit),'k')
% legend('entrada A/D','Saida FPGA')
% title('treinamento imaginario')
% 
% figure 
% hold on
% plot(0:2/length(S_DD_Fit):2-2/length(S_DD_Fit),db(abs(fft(S_DD_Fit))),'r')
% plot(0:2/length(S_A_Fit):2-2/length(S_A_Fit),db(abs(fft(S_A_Fit))),'k')
% legend('entrada A/D','Saida FPGA')
% title('treinamento')
% 
% figure 
% hold on
% plot(S_DD_Fiv,'r')
% plot(S_A_Fiv,'k')
% legend('entrada A/D','Saida FPGA')
% title('validacao')
% 
% figure 
% hold on
% plot(0:2/length(S_DD_Fit):2-2/length(S_DD_Fit),db(abs(fft(S_DD_Fiv))),'r')
% plot(0:2/length(S_A_Fit):2-2/length(S_A_Fit),db(abs(fft(S_A_Fiv))),'k')
% legend('entrada A/D','Saida FPGA')
% title('validacao')
% 
% figure 
% hold on
% plot(S_DD_Fis,'r')
% plot(S_A_Fis,'k')
% legend('entrada A/D','Saida FPGA')
% title('sinal')
% 
% figure 
% hold on
% plot(0:2/length(S_DD_Fit):2-2/length(S_DD_Fit),db(abs(fft(S_DD_Fis))),'r')
% plot(0:2/length(S_A_Fit):2-2/length(S_A_Fit),db(abs(fft(S_A_Fis))),'k')
% legend('entrada A/D','Saida FPGA')
% title('sinal')

% xlim([100 110])
max_m2=(max_delay-1)*2+1;%impar


S_DD_Fita=S_DD_Fit;
S_A_Fita=S_A_Fit;
S_DD_Fiva=S_DD_Fiv;
S_A_Fiva=S_A_Fiv;
S_DD_Fisa=S_DD_Fis;
S_A_Fisa=S_A_Fis;
Tam_a=length(S_DD_Fita);

if procura_atr_fino==1 
lim_S_m2=max_m2;
else 
lim_S_m2=1;
end

for m2=1:1:lim_S_m2


m2_c=m2-(max_m2-1)/2-1;

N_P=length(S_DD_Fit);


% figure
% hold on
% plot(real(S_A_Fiv(ini_nmse:end)),'k')
% plot(real(S_DD_Fiv(ini_nmse:end)),'r')


if procura_atr_fino==1 
lim_I_delay_a=ini_t_delay;
lim_S_delay_a=max_delay;
else 
lim_I_delay_a=m2_c+1;
lim_S_delay_a=m2_c+1;
end

for f_fir=1:1:1

    
for ordem=ini_u_ordem:1:max_u_ordem
for delay_a=ini_t_delay:passo_delay:max_delay
for KarnelV=iKarnelV:1:min([fKarnelV ordem delay_a])
    
m=delay_a-1; 

if procura_atr_fino==1 
% %%%%%%%%%%%%%%%%%%%%%%%%% mudar aqui
Sx=S_A_Fit(1+m+m2_c+max_delay:-max_delay+N_P);
Sy=S_DD_Fit(1+max_delay:-max_delay+N_P-m-m2_c);
Sx_v=S_A_Fiv(1+m+m2_c+max_delay:-max_delay+N_P);
Sy_v=S_DD_Fiv(1+max_delay:-max_delay+N_P-m-m2_c);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% aqui
else  
Sx=S_A_Fit;
Sy=S_DD_Fit;
Sx_v=S_A_Fiv;
Sy_v=S_DD_Fiv;
end
    
% figure
% hold on
% plot(real(Sx))
% plot(real(Sy),'r')

ini=1;

Ux=[];
Uy=[];
N_PP=length(Sx);

for k=0:m
    Ux(:,k+1)=Sx(m+1-k:N_PP-k,1);
end
for k=1:m
    Uy(:,k)=Sy(m+1-k:N_PP-k,1);
end

U=[Ux(:,ini:end)];
for k=1:1:ordem-1
U=[U Ux(:,ini:end).*abs(Ux(:,ini:end)).^(k*k_t_amp_pot)];
end

if sem_DC==0
U=[U ones(size(U(:,1)))];
else
U=[U];
end

if sem_REALIM==0
U=[U Uy];
else
U=[U];
end

if KarnelV>=2
for kKarnelV=2:1:KarnelV
[VetComR,VetSemR,NVetComR,NVetSemR]=Gera_coef_Volterra(kKarnelV);
Ua=[];
for k9=1:1:NVetSemR(1)
Sa=ones(length(U),1);
for k8=1:1:NVetSemR(2)
Sa=Sa.*Sx(m+1-(VetSemR(k9,k8)-1):N_PP-(VetSemR(k9,k8)-1),1);  
end
Ua=[Ua Sa];
end
U=[U Ua];
end
end

Meu_y=Sy(1+m:N_PP);

Mat_Acmp=[U Meu_y];
% Theta=inv(U'*U)*U'*Meu_y;

figure
hold on
plot(real(Meu_y))

Theta=U\Meu_y;

% figure
% hold on
% plot(S_A_Fit(ini_nmse:end),'m')
% plot(S_DD_Fit(ini_nmse:end),'c')
% legend('A','DD')
% title(num2str(Theta'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55555 Continuar daqui
for k=1:1:length(Theta)
    if abs(Theta(k))<1*10^-lim_casas
        Theta(k)=0;
    end
end

ae2_1=zeros(1,delay_a);
ae2_2=zeros(1,delay_a);
ae2_3=zeros(1,delay_a);
ae2_4=zeros(1,delay_a);
ae2_5=zeros(1,delay_a);
ae2_6=zeros(1,delay_a);
ae2_7=zeros(1,delay_a);
ae2_8=zeros(1,delay_a);

if ordem>0
ae2_1=transpose(Theta(1:delay_a));
end
if ordem>1
ae2_2=transpose(Theta(delay_a+1:2*delay_a));
end
if ordem>2
ae2_3=transpose(Theta(2*delay_a+1:3*delay_a));
end
if ordem>3
ae2_4=transpose(Theta(3*delay_a+1:4*delay_a));
end
if ordem>4
ae2_5=transpose(Theta(4*delay_a+1:5*delay_a));
end
if ordem>5
ae2_6=transpose(Theta(5*delay_a+1:6*delay_a));
end
if ordem>6
ae2_7=transpose(Theta(6*delay_a+1:7*delay_a));
end
if ordem>7
ae2_8=transpose(Theta(7*delay_a+1:8*delay_a));
end


if sem_REALIM==0
be2_1=[1 -transpose(Theta(end-delay_a+2:end))];
else
be2_1=[1,zeros(1,delay_a-1)];
ao=Theta(ordem*delay_a+1:end,1);
end

if sem_DC==0
libe2_1=Theta(end-delay_a+1)/sum(be2_1);
else
libe2_1=0;    
end
% disp(['a_n=          ',num2str(ae2_1),''])
% disp(['b_n=          ',num2str(be2_1),''])


be4_1=be2_1;
ae4_1=ae2_1/(be4_1(1));
ae4_2=ae2_2/(be4_1(1));
ae4_3=ae2_3/(be4_1(1));
ae4_4=ae2_4/(be4_1(1));
ae4_5=ae2_5/(be4_1(1));
ae4_6=ae2_6/(be4_1(1));
ae4_7=ae2_7/(be4_1(1));
ae4_8=ae2_8/(be4_1(1));

libe4_1=libe2_1;

for k_o=1:1:length(ae4_1)
    a_m1(ordem,delay_a,KarnelV,m2,k_o)=ae4_1(k_o);
    a_m2(ordem,delay_a,KarnelV,m2,k_o)=ae4_2(k_o);
    a_m3(ordem,delay_a,KarnelV,m2,k_o)=ae4_3(k_o);
    a_m4(ordem,delay_a,KarnelV,m2,k_o)=ae4_4(k_o);
    a_m5(ordem,delay_a,KarnelV,m2,k_o)=ae4_5(k_o);
    a_m6(ordem,delay_a,KarnelV,m2,k_o)=ae4_6(k_o);
    a_m7(ordem,delay_a,KarnelV,m2,k_o)=ae4_7(k_o);
    a_m8(ordem,delay_a,KarnelV,m2,k_o)=ae4_8(k_o);
end
for k_o=1:1:length(be4_1)
    b_m(ordem,delay_a,KarnelV,m2,k_o)=be4_1(k_o);
end
for k_o=1:1:length(ao)
    ao_m(ordem,delay_a,KarnelV,m2,k_o)=ao(k_o);
end
sz_ao(ordem,delay_a,KarnelV,m2)=length(ao);
KarnelV_m(ordem,delay_a,KarnelV,m2)=KarnelV;
lib_m(ordem,delay_a,KarnelV,m2)=libe4_1;
Sy4=[];


% figure
% hold on
% plot(real(Sx_v))

Sx4=Sx_v;
for k=max([length(ae4_1) length(be4_1)]):1:length(Sx4)
   Sy4(k)=0;
    for ka=1:1:length(ae4_1)
        Sy4(k)=Sy4(k)+Sx4(k-ka+1)*ae4_1(ka)+Sx4(k-ka+1)*abs(Sx4(k-ka+1)).^(1*k_t_amp_pot)*ae4_2(ka)+Sx4(k-ka+1)*abs(Sx4(k-ka+1)).^(2*k_t_amp_pot)*ae4_3(ka)...
            +Sx4(k-ka+1)*abs(Sx4(k-ka+1)).^(3*k_t_amp_pot)*ae4_4(ka)+Sx4(k-ka+1)*abs(Sx4(k-ka+1)).^(4*k_t_amp_pot)*ae4_5(ka)+Sx4(k-ka+1)*abs(Sx4(k-ka+1)).^(5*k_t_amp_pot)*ae4_6(ka)...
            +Sx4(k-ka+1)*abs(Sx4(k-ka+1)).^(6*k_t_amp_pot)*ae4_7(ka)+Sx4(k-ka+1)*abs(Sx4(k-ka+1)).^(7*k_t_amp_pot)*ae4_8(ka);
    end
    for kb=2:1:length(be4_1) %começa do segundo termo 
        Sy4(k)=Sy4(k)-Sy4(k-kb+1)*be4_1(kb);
    end
end
Sy4=Sy4+libe4_1;


kKar=0;
for kKarnelV=2:1:KarnelV
Sy4(kKarnelV,:)=zeros(1,length(Sx4));
[VetComR,VetSemR,NVetComR,NVetSemR]=Gera_coef_Volterra(kKarnelV);


for k=1:1:NVetSemR(1)
   kKar=kKar+1;
   a(k)=ao(kKar); 
end

for k=kKarnelV:1:length(Sx4)
    for k2=1:1:NVetSemR(1)
        Sx4A=1;
        for k3=1:1:kKarnelV
            Sx4A=Sx4A*Sx4(k+1-VetSemR(k2,k3));
        end
        Sy4(kKarnelV,k)=Sy4(kKarnelV,k)+a(k2)*Sx4A;
    end
end
end

[tam1,tam2]=size(Sy4);
if tam1>1
    Sy4=sum(Sy4);
end


[nmse_dB(ordem,delay_a,KarnelV,m2)]=test_nmse(transpose(Sy_v(ini_nmse:end)), Sy4(ini_nmse:end));

MED1(ordem,delay_a,KarnelV,m2)=gfit2(transpose(Sy_v(ini_nmse:end)), Sy4(ini_nmse:end),'1');
MED2(ordem,delay_a,KarnelV,m2)=gfit2(transpose(Sy_v(ini_nmse:end)), Sy4(ini_nmse:end),'2');
MED3(ordem,delay_a,KarnelV,m2)=gfit2(transpose(Sy_v(ini_nmse:end)), Sy4(ini_nmse:end),'3');
MED4(ordem,delay_a,KarnelV,m2)=gfit2(transpose(Sy_v(ini_nmse:end)), Sy4(ini_nmse:end),'4');
MED5(ordem,delay_a,KarnelV,m2)=gfit2(transpose(Sy_v(ini_nmse:end)), Sy4(ini_nmse:end),'5');
MED6(ordem,delay_a,KarnelV,m2)=gfit2(transpose(Sy_v(ini_nmse:end)), Sy4(ini_nmse:end),'6');
MED7(ordem,delay_a,KarnelV,m2)=gfit2(transpose(Sy_v(ini_nmse:end)), Sy4(ini_nmse:end),'7');
MED8(ordem,delay_a,KarnelV,m2)=gfit2(transpose(Sy_v(ini_nmse:end)), Sy4(ini_nmse:end),'8');
MED9(ordem,delay_a,KarnelV,m2)=gfit2(transpose(Sy_v(ini_nmse:end)), Sy4(ini_nmse:end),'9');
MED10(ordem,delay_a,KarnelV,m2)=gfit2(transpose(Sy_v(ini_nmse:end)), Sy4(ini_nmse:end),'10');
MED11(ordem,delay_a,KarnelV,m2)=gfit2(transpose(Sy_v(ini_nmse:end)), Sy4(ini_nmse:end),'11');

disp(['delaya=',num2str(delay_a),' ordem=',num2str(ordem),' KarnelV=',num2str(KarnelV),' m2=',num2str(m2),' max_m2=',num2str(max_m2)])
end
end
end
end
end

nmse_dB;
% figure
% hold on

nmse_dB3=[];
MED1_3=[];
MED2_3=[];
MED3_3=[];
MED4_3=[];
MED5_3=[];
MED6_3=[];
MED7_3=[];
MED8_3=[];
MED9_3=[];
MED10_3=[];
MED11_3=[];

min_nmse=9999999;
delay_a3=999999;
nmse_M=[];

for k3=iKarnelV:1:fKarnelV
for d3=ini_t_delay:1:max_delay
for o3=ini_u_ordem:1:max_u_ordem
for m3=1:1:lim_S_m2
        nmse_dB3(m3)=nmse_dB(o3,d3,k3,m3);
        MED1_3(m3)=MED1(o3,d3,k3,m3);
        MED2_3(m3)=MED2(o3,d3,k3,m3);
        MED3_3(m3)=MED3(o3,d3,k3,m3);
        MED4_3(m3)=MED4(o3,d3,k3,m3);
        MED5_3(m3)=MED5(o3,d3,k3,m3);
        MED6_3(m3)=MED6(o3,d3,k3,m3);
        MED7_3(m3)=MED7(o3,d3,k3,m3);
        MED8_3(m3)=MED8(o3,d3,k3,m3);
        MED9_3(m3)=MED9(o3,d3,k3,m3);
        MED10_3(m3)=MED10(o3,d3,k3,m3);
        MED11_3(m3)=MED11(o3,d3,k3,m3);
        D_3(m3)=d3;
        if nmse_dB3(m3)<min_nmse
            if nmse_dB3(m3)<-150 && D_3(m3)<delay_a3
                min_nmse=nmse_dB3(m3);
                delay_a3=D_3(m3);
                m_m3=m3;
                o_o3=o3;
            else
                if nmse_dB3(m3)<min_nmse
                    min_nmse=nmse_dB3(m3);
                    delay_a3=D_3(m3);
                    m_m3=m3;
                    o_o3=o3;
                    k_k3=k3;
                end
            end
        end
end
[delay_M(o3,d3,k3),posB]=min(nmse_dB3);
delay_MED1(o3,d3,k3)=MED1_3(posB);
delay_MED2(o3,d3,k3)=MED2_3(posB);
delay_MED3(o3,d3,k3)=MED3_3(posB);
delay_MED4(o3,d3,k3)=MED4_3(posB);
delay_MED5(o3,d3,k3)=MED5_3(posB);
delay_MED6(o3,d3,k3)=MED6_3(posB);
delay_MED7(o3,d3,k3)=MED7_3(posB);
delay_MED8(o3,d3,k3)=MED8_3(posB);
delay_MED9(o3,d3,k3)=MED9_3(posB);
delay_MED10(o3,d3,k3)=MED10_3(posB);
delay_MED11(o3,d3,k3)=MED11_3(posB);
end
end
end

% for d4=1:1:max_delay+1
% for o4=1:1:max_u_ordem+1
% 
% if o4==1 && d4==1
% delay_MN(o4,d4)=0;
% end  
% if o4==1 && d4>1
% delay_MN(o4,d4)=d4-1;
% end  
% if o4>1 && d4==1
% delay_MN(o4,d4)=o4-1;
% end  
% if o4>1 && d4>1
% delay_MN(o4,d4)=delay_M(o4-1,d4-1);
% end
% 
% end
% end

disp(delay_M)
disp('mse')
disp(delay_MED1)
disp('nmse')
disp(delay_MED2)
disp('rmse')
disp(delay_MED3)
disp('nrmse')
disp(delay_MED4)
disp('mae')
disp(delay_MED5)
disp('mare')
disp(delay_MED6)
disp('coefficient of correlation (r)')
disp(delay_MED7)
disp('coefficient of determination (d)')
disp(delay_MED8)
disp('coefficient of efficiency (e)')
disp(delay_MED9)
disp('maximum absolute error')
disp(delay_MED10)
disp('maximum absolute relative error')
disp(delay_MED11)
% plot(nmse_dB(o_o3,delay_a3,:),'k')
% plot(nmse_dB(2,:),'b')
% Ylim([-160 0])

% [y_fir,x_fir]=min(nmse_dB(1,:));
% [y_iir,x_iir]=min(nmse_dB(2,:));

for D_i=1:delay_a3
    a_fir1(D_i)=a_m1(o_o3,delay_a3,k_k3,m_m3,D_i);
    a_fir2(D_i)=a_m2(o_o3,delay_a3,k_k3,m_m3,D_i);
    a_fir3(D_i)=a_m3(o_o3,delay_a3,k_k3,m_m3,D_i);
    a_fir4(D_i)=a_m4(o_o3,delay_a3,k_k3,m_m3,D_i);
    a_fir5(D_i)=a_m5(o_o3,delay_a3,k_k3,m_m3,D_i);
    a_fir6(D_i)=a_m6(o_o3,delay_a3,k_k3,m_m3,D_i);
    a_fir7(D_i)=a_m7(o_o3,delay_a3,k_k3,m_m3,D_i);
    a_fir8(D_i)=a_m8(o_o3,delay_a3,k_k3,m_m3,D_i);
    b_fir(D_i)=b_m(o_o3,delay_a3,k_k3,m_m3,D_i);
end

for D_i=1:1:sz_ao(o_o3,delay_a3,k_k3,m_m3)
    ao_fir(D_i)=ao_m(o_o3,delay_a3,k_k3,m_m3,D_i);
end

KarnelV_fir=KarnelV_m(o_o3,delay_a3,k_k3,m_m3);
lib_fir=lib_m(o_o3,delay_a3,k_k3,m_m3);

% a_iir(1,1:x_iir)=a_m(2,x_iir,1:x_iir);
% b_iir(1,1:x_iir)=b_m(2,x_iir,1:x_iir);
% lib_iir=lib_m(2,x_iir);
m=delay_a3-1;
m3_c=m_m3-(max_m2-1)/2-1;

% figure 
% hold on
% plot(S_DD_Fisa,'r')
% plot(S_A_Fisa,'k')
% legend('Depois','Antes')
% title('Antes')

% S_DD_Fis=S_DD_Fisa(1+max_m2+m3_c:end-max_m2+m3_c-2);
% S_A_Fis=S_A_Fisa(1+max_m2:end-max_m2);
% 
% figure 
% hold on
% plot(S_DD_Fis,'r')
% plot(S_A_Fis,'k')
% title('Ajustado')

% Sx_v=S_A_Fiv(1+m+m2_c+max_delay:-max_delay+N_P);
% Sy_v=S_DD_Fiv(1+max_delay:-max_delay+N_P-m-m2_c);
if procura_atr_fino==1 
Sx_s=S_A_Fis(1+m+m3_c+max_delay:-max_delay+N_P);
Sy_s=S_DD_Fis(1+max_delay:-max_delay+N_P-m-m3_c);
else
Sx_s=S_A_Fis;
Sy_s=S_DD_Fis;
end
Sx_fir=Sx_s;
for k=max([length(a_fir1) length(b_fir)]):1:length(Sx_fir)
   Sy_fir(k)=0;
    for ka=1:1:length(a_fir1)
        Sy_fir(k)=Sy_fir(k)+Sx_fir(k-ka+1)*a_fir1(ka)+Sx_fir(k-ka+1)*abs(Sx_fir(k-ka+1)).^(1*k_t_amp_pot)*a_fir2(ka)+Sx_fir(k-ka+1)*abs(Sx_fir(k-ka+1)).^(2*k_t_amp_pot)*a_fir3(ka)...
            +Sx_fir(k-ka+1)*abs(Sx_fir(k-ka+1)).^(3*k_t_amp_pot)*a_fir4(ka)+Sx_fir(k-ka+1)*abs(Sx_fir(k-ka+1)).^(4*k_t_amp_pot)*a_fir5(ka)...
            +Sx_fir(k-ka+1)*abs(Sx_fir(k-ka+1)).^(5*k_t_amp_pot)*a_fir6(ka)+Sx_fir(k-ka+1)*abs(Sx_fir(k-ka+1)).^(6*k_t_amp_pot)*a_fir7(ka)...
            +Sx_fir(k-ka+1)*abs(Sx_fir(k-ka+1)).^(7*k_t_amp_pot)*a_fir8(ka);
    end
    for kb=2:1:length(b_fir) %começa do segundo termo 
        Sy_fir(k)=Sy_fir(k)-Sy_fir(k-kb+1)*b_fir(kb);
    end
end
Sy_fir=Sy_fir+lib_fir;

kKar=0;
for kKarnelV=2:1:KarnelV_fir
Sy_fir(kKarnelV,:)=zeros(1,length(Sx_fir));
[VetComR,VetSemR,NVetComR,NVetSemR]=Gera_coef_Volterra(kKarnelV);


for k=1:1:NVetSemR(1)
   kKar=kKar+1;
   a(k)=ao_fir(kKar); 
end

for k=kKarnelV:1:length(Sx_fir)
    for k2=1:1:NVetSemR(1)
        Sx4A=1;
        for k3=1:1:kKarnelV
            Sx4A=Sx4A*Sx_fir(k+1-VetSemR(k2,k3));
        end
        Sy_fir(kKarnelV,k)=Sy_fir(kKarnelV,k)+a(k2)*Sx4A;
    end
end
end

[tam1,tam2]=size(Sy_fir);
if tam1>1
    Sy_fir=sum(Sy_fir);
end

% Sx_iir=Sx_s;
% for k=max([length(a_iir) length(b_iir)]):1:length(Sx_iir)
%    Sy_iir(k)=0;
%     for ka=1:1:length(a_iir)
%         Sy_iir(k)=Sy_iir(k)+Sx_iir(k-ka+1)*a_iir(ka);
%     end
%     for kb=2:1:length(b_iir) %começa do segundo termo 
%         Sy_iir(k)=Sy_iir(k)-Sy_iir(k-kb+1)*b_iir(kb);
%     end
% end
% Sy_iir=Sy_iir+lib_iir;

if procura_atr_fino~=0
    Sx_s=S_A_Fis(1+max_delay:-max_delay+N_P-m-m3_c);
end
% figure
% hold on
% plot(Sx_s(ini_nmse:end),'c')
% plot(Sy_s(ini_nmse:end),'r')
% plot(Sy_fir(ini_nmse:end),'k')
% % plot(Sy_iir(ini_nmse:end),'b')
% legend('Sinal_D','certo','modelo fir')%,'modelo iir')

[nmse_fir]=test_nmse(transpose(Sy_s(ini_nmse:end)),Sy_fir(ini_nmse:end));
% [nmse_iir]=test_nmse(transpose(Sy_s(ini_nmse:end)),Sy_iir(ini_nmse:end));
[nmse_SF]=test_nmse(Sy_s(ini_nmse:end),Sx_s(ini_nmse:end));

% figure
% hold on
% plot(0:2/length(Sx_s(ini_nmse:end)):2-2/length(Sx_s(ini_nmse:end)),db(abs(fft(Sx_s(ini_nmse:end)))),'c')
% plot(0:2/length(Sy_s(ini_nmse:end)):2-2/length(Sy_s(ini_nmse:end)),db(abs(fft(Sy_s(ini_nmse:end)))),'r')
% plot(0:2/length(Sy_fir(ini_nmse:end)):2-2/length(Sy_fir(ini_nmse:end)),db(abs(fft(Sy_fir(ini_nmse:end)))),'k')
% % plot(Sy_iir(ini_nmse:end),'b')
% legend('Sinal_D','certo','modelo fir')%,'modelo iir')


disp(['nmse_fir =', num2str(nmse_fir),' nmse_SF =', num2str(nmse_SF)])
disp('    ')

if bandaS2(1)==0 && bandaS2(2)==1
A=1;
B=1;
elseif bandaS2(1)==0 && bandaS2(2)<1
[B,A] = BUTTER(ordem_C2,bandaS2(2),'low');    
elseif bandaS2(1)>0 && bandaS2(2)==1
[B,A] = BUTTER(ordem_C2,bandaS2(1),'high');    
elseif bandaS2(1)>0 && bandaS2(2)<1
[B,A] = BUTTER(ordem_C2,bandaS2);  
end

Sx_s = filter(B,A,Sx_s);
Sy_s = filter(B,A,Sy_s);
Sy_fir = filter(B,A,Sy_fir);

% figure
% hold on
% plot(real(Sx_s(ini_nmse:end)),'c')
% plot(real(Sy_s(ini_nmse:end)),'r')
% plot(real(Sy_fir(ini_nmse:end)),'k')
% % plot(Sy_iir(ini_nmse:end),'b')
% legend('Sinal_D','certo','modelo fir')%,'modelo iir')
% 
% 
% figure
% hold on
% plot(real(Sx_s(ini_nmse:end)),'c')
% plot(real(Sy_s(ini_nmse:end)),'r')
% plot(real(Sy_fir(ini_nmse:end)),'k')
% % plot(Sy_iir(ini_nmse:end),'b')
% legend('Sinal_D','certo','modelo fir')%,'modelo iir')
% 


[nmse_fir]=test_nmse(transpose(Sy_s(ini_nmse:end)),Sy_fir(ini_nmse:end));

gfit2(transpose(Sy_s(ini_nmse:end)),Sy_fir(ini_nmse:end),'1')
gfit2(transpose(Sy_s(ini_nmse:end)),Sy_fir(ini_nmse:end),'2')
gfit2(transpose(Sy_s(ini_nmse:end)),Sy_fir(ini_nmse:end),'3')
gfit2(transpose(Sy_s(ini_nmse:end)),Sy_fir(ini_nmse:end),'4')
gfit2(transpose(Sy_s(ini_nmse:end)),Sy_fir(ini_nmse:end),'5')
gfit2(transpose(Sy_s(ini_nmse:end)),Sy_fir(ini_nmse:end),'6')
gfit2(transpose(Sy_s(ini_nmse:end)),Sy_fir(ini_nmse:end),'7')
gfit2(transpose(Sy_s(ini_nmse:end)),Sy_fir(ini_nmse:end),'8')
gfit2(transpose(Sy_s(ini_nmse:end)),Sy_fir(ini_nmse:end),'9')
gfit2(transpose(Sy_s(ini_nmse:end)),Sy_fir(ini_nmse:end),'10')
gfit2(transpose(Sy_s(ini_nmse:end)),Sy_fir(ini_nmse:end),'11')

% [nmse_iir]=test_nmse(transpose(Sy_s(ini_nmse:end)),Sy_iir(ini_nmse:end));
[nmse_SF]=test_nmse(Sy_s(ini_nmse:end),Sx_s(ini_nmse:end));

% figure
% hold on
% plot(0:2/length(Sx_s(ini_nmse:end)):2-2/length(Sx_s(ini_nmse:end)),db(abs(fft(Sx_s(ini_nmse:end)))),'c')
% plot(0:2/length(Sy_s(ini_nmse:end)):2-2/length(Sy_s(ini_nmse:end)),db(abs(fft(Sy_s(ini_nmse:end)))),'r')
% plot(0:2/length(Sy_fir(ini_nmse:end)):2-2/length(Sy_fir(ini_nmse:end)),db(abs(fft(Sy_fir(ini_nmse:end)))),'k')
% % plot(Sy_iir(ini_nmse:end),'b')
% legend('Sinal_D','certo','modelo fir')%,'modelo iir')


Hs = spectrum.welch('Hamming',7000);
Fs=76.8*10^6;

% figure
% hold on
% subplot(1,3,1)
% psd(Hs,Sx_s,'Fs',Fs,'CenterDC',true)
% legend('Sinal_D')%,'modelo iir')
% ylim([-160 -70])
% xlim([-8 8])
% subplot(1,3,2)
% psd(Hs,Sy_s,'Fs',Fs,'CenterDC',true)
% legend('certo')%,'modelo iir')
% ylim([-160 -70])
% xlim([-8 8])
% subplot(1,3,3)
% psd(Hs,Sy_fir,'Fs',Fs,'CenterDC',true)
% legend('modelo fir')%,'modelo iir')
% ylim([-160 -70])
% xlim([-8 8])

% figure
% hold on
% psd(Hs,Sx_s,'Fs',Fs,'CenterDC',true)
% psd(Hs,Sy_s,'Fs',Fs,'CenterDC',true)
% psd(Hs,Sy_fir,'Fs',Fs,'CenterDC',true)


disp(['nmse_modelo =', num2str(nmse_fir),' nmse_SF =', num2str(nmse_SF)])
disp('    ')
if o_o3>=1
disp(['a1 = ',num2str((a_fir1)*k_mul_a^0)])
end
if o_o3>=2
disp(['a2 = ',num2str((a_fir2)*k_mul_a^0)])
end
if o_o3>=3
disp(['a3 = ',num2str((a_fir3)*k_mul_a^0)])
end
if o_o3>=4
disp(['a4 = ',num2str((a_fir4)*k_mul_a^0)])
end
if o_o3>=5
disp(['a5 = ',num2str((a_fir5)*k_mul_a^0)])
end
if o_o3>=6
disp(['a6 = ',num2str((a_fir6)*k_mul_a^0)])
end
if o_o3>=7
disp(['a7 = ',num2str((a_fir7)*k_mul_a^0)])
end
if o_o3>=8
disp(['a8 = ',num2str((a_fir8)*k_mul_a^0)])
end

if KarnelV_fir>=2
disp(['ao2 = ',num2str((ao_fir(1:1))*k_mul_a^0)])
end
if KarnelV_fir>=3
disp(['ao3 = ',num2str((ao_fir(1+1:1+7))*k_mul_a^0)])
end
if KarnelV_fir>=4
disp(['ao4 = ',num2str((ao_fir(1+7+1:1+7+31))*k_mul_a^0)])
end
if KarnelV_fir>=5
disp(['ao5 = ',num2str((ao_fir(1+7+31+1:1+7+31+121))*k_mul_a^0)])
end
if KarnelV_fir>=6
disp(['ao6 = ',num2str((ao_fir(1+7+31+121+1:1+7+31+121+456))*k_mul_a^0)])
end
if KarnelV_fir>=7
disp(['ao7 = ',num2str((ao_fir(1+7+31+121+456+1:1+7+31+121+456+1709))*k_mul_a^0)])
end
if KarnelV_fir>=8
disp(['ao8 = ',num2str((ao_fir(1+7+31+121+456+1709+1:1+7+31+121+456+1709+6427))*k_mul_a^0)])
end

% disp(['a1 = ',num2str(real(a_fir1)*k_mul_a^0)])
% disp(['a1i= ',num2str(imag(a_fir1)*k_mul_a^0)])
% disp(['a2 = ',num2str(real(a_fir2)*k_mul_a^0)])
% disp(['a2i= ',num2str(imag(a_fir2)*k_mul_a^0)])
% disp(['a3 = ',num2str(real(a_fir3)*k_mul_a^0)])
% disp(['a3i= ',num2str(imag(a_fir3)*k_mul_a^0)])
% disp(['a4 = ',num2str(real(a_fir4)*k_mul_a^0)])
% disp(['a4i= ',num2str(imag(a_fir4)*k_mul_a^0)])
% disp(['a5 = ',num2str(real(a_fir5)*k_mul_a^0)])
% disp(['a5i= ',num2str(imag(a_fir5)*k_mul_a^0)])
% disp(['a6 = ',num2str(real(a_fir6)*k_mul_a^0)])
% disp(['a6i= ',num2str(imag(a_fir6)*k_mul_a^0)])
% disp(['a7 = ',num2str(real(a_fir7)*k_mul_a^0)])
% disp(['a7i= ',num2str(imag(a_fir7)*k_mul_a^0)])
% disp(['a8 = ',num2str(real(a_fir8)*k_mul_a^0)])
% disp(['a8i= ',num2str(imag(a_fir8)*k_mul_a^0)])
% 
% disp(['b = ',num2str(real(b_fir)*k_mul_a^0)])
% disp(['bi= ',num2str(imag(b_fir)*k_mul_a^0)])
% disp(['lib= ',num2str(lib_fir*k_mul_a^0)])

Fs=1.25*2*76.8*10^6;

% 
% figF2=figure;
% hold on
% grid on
% psdSx_s=psd(Hs,Sx_s,'Fs',Fs,'CenterDC',true);
% psdSy_s=psd(Hs,Sy_s,'Fs',Fs,'CenterDC',true);
% psdSy_fir=psd(Hs,Sy_fir,'Fs',Fs,'CenterDC',true);
% plot((get(psdSx_s,'Frequencies')),db(get(psdSx_s,'Data'))/2,'c')
% plot((get(psdSy_s,'Frequencies')),db(get(psdSy_s,'Data'))/2,'r')
% plot((get(psdSy_fir,'Frequencies')),db(get(psdSy_fir,'Data'))/2,'k')
% ylim([-160 -100])
% xlim([-21 21]*010^6)
% legend('Sinal sem amplificador','Sinal com amplificador','modelo')
% title(['NMSE sem modelo',num2str(nmse_SF),' Modelo NMSE:',num2str(nmse_fir),' o:',num2str(o_o3),' d:',num2str(length(a_fir1))])


figF2=figure;
hold on
% grid on
AmpFIN=15000;
Sx_s=Sx_s*AmpFIN;Sy_s=Sy_s*AmpFIN;Sy_fir=Sy_fir*AmpFIN;
psdSx_s=psd(Hs,Sx_s,'Fs',Fs,'CenterDC',true);
psdSy_s=psd(Hs,Sy_s,'Fs',Fs,'CenterDC',true);
psdSy_fir=psd(Hs,Sy_fir,'Fs',Fs,'CenterDC',true);
% plot((get(psdSx_s,'Frequencies')),db(get(psdSx_s,'Data'))/2,'color',[0,0,0],'LineWidth',2)
plot((get(psdSy_s,'Frequencies')),db(get(psdSy_s,'Data'))/2,'color',[0,0,0],'LineWidth',2)
plot((get(psdSy_fir,'Frequencies')),db(get(psdSy_fir,'Data'))/2,'color',[0.4,0.4,0.4],'LineWidth',2)
ylim([-80 0])
xlim([-27 27]*010^6)
set(gca,'FontSize',15)
ylabel('Amplitude (dB)')
xlabel('Frequência')
legend('Sinsl com Amplificador','Modelo','Location','South')
% title(['NMSE (Sinal recebido x Modelo):',num2str(nmse_fir),' NMSE (Sinal recebido x Sinal Enviado):',num2str(nmse_SF),' order:',num2str(o_o3),' delay:',num2str(length(a_fir1))])
% title('Sinal com ordem 7 e 4 de atraso')
title([num2str(k_Norm_k),' ',num2str(k_max_o_k),' ',num2str(k_t_amp_pot),' nmse_modelo =', num2str(nmse_fir),' nmse_SF =', num2str(nmse_SF)])

vet_nmse_t_amp_pot_k=[vet_nmse_t_amp_pot_k (nmse_fir)];
vet_nmse_t_amp_pot_o=[vet_nmse_t_amp_pot_o k_Norm_k];



end
end



plot(vet_nmse_t_amp_pot_o,vet_nmse_t_amp_pot_k)
title(['Norm ',num2str(k_Norm_k),'O ',num2str(k_max_o_k),'D ',num2str(k_max_d_k),'amp ',num2str(k_t_amp_pot)])

end
end
end
end