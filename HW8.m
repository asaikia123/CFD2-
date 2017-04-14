clear all;
clc;
format long
e=1;
de=0.01837;
M=2;
V=1;
rho=1;
u=V*cosd(15);

dn=tand(40)/(26-1);
n= 0 :dn:tand(40);
gamma=1.4;

E=zeros(2,length(n),3);
F=zeros(2,length(n),3);

v=0;
p=rho*((1/(gamma*M^2))+ ((gamma-1)/(2*gamma))*(1-u^2-v^2));

%% at wall %%
E(1,1,1)=rho*u;
E(1,1,2)=rho*u^2+p;
E(1,1,3)=rho*u*v;
F(1,1,1)=rho*v;
F(1,1,3)=rho*v^2+p;
F(1,1,2)=rho*u*v; 



v=-V*sind(15);

p=rho*((1/(gamma*M^2))+ ((gamma-1)/(2*gamma))*(1-u^2-v^2));


E(1,2,1)=rho*u;
E(1,2,2)=rho*u^2+p;
E(1,2,3)=rho*u*v;
F(1,2,1)=rho*v;
F(1,2,3)=rho*v^2+p;
F(1,2,2)=rho*u*v;
Ht=(gamma/(gamma-1)) *( E(1,1,2)-(E(1,1,1)*u)/(E(1,1,1)/u)) +(u^2+v^2)/2 ;

E_bar=zeros(2,length(n),3);
F_bar=zeros(2,length(n),3);


%% Encoding %%
%% MacCormac Scheme %% 
i=1;
for j=2:length(n)-1;


    %% Predictor %%
E_bar(i+1,j,:)=E(i,j,:)-(de/dn)*(F(i,j+1,:)-F(i,j,:));
%% decoding for predictor %%
v=E_bar(i+1,j,3)/E_bar(i+1,j,1);
Ht=(gamma/(gamma-1)) *( E_bar(i+1,j,2)-(E_bar(i+1,j,1)*u)/(E_bar(i+1,j,1)/u)) +(u^2+v^2)/2;
u=(gamma/(gamma+1))*(E_bar(i+1,j,2)/E_bar(i+1,j,1))+  sqrt( ((gamma*E_bar(i+1,j,2))/((gamma+1)*E_bar(i+1,j,1)))^2 - ((gamma-1)/(gamma+1)) *(2*Ht-v^2));
p=E_bar(i+1,j,2)-u*E_bar(i+1,j,1);
rho=E_bar(i+1,j,1)/u;

%% encoding for predictor %%
F_bar(i+1,j,1)=rho*v;
F_bar(i+1,j,3)=rho*v^2+p;
F_bar(i+1,j,2)=rho*u*v;

%%  Corrector

 E(i+1,j,:)=0.5*(( E(i,j,:)+E_bar(i+1,j,:))   - (de/dn)*(F_bar(i+1,j,:)-F_bar(i+1,j-1,:)));
 
 %% decoding for predictor %%
v=E(i+1,j,3)/E(i+1,j,1);
Ht=(gamma/(gamma-1)) *( E(i+1,j,2)-(E(i+1,j,1)*u)/(E(i+1,j,1)/u)) +(u^2+v^2)/2;
u=(gamma/(gamma+1))*(E(i+1,j,2)/E(i+1,j,1))+  sqrt( ((gamma*E_bar(i+1,j,2))/((gamma+1)*E(i+1,j,1)))^2 - ((gamma-1)/(gamma+1)) *(2*Ht-v^2));
p=E(i+1,j,2)-u*E(i+1,j,1);
rho=E(i+1,j,1)/u;

%% encoding for predictor %%
F(i+1,j,1)=rho*v;
F(i+1,j,3)=rho*v^2+p;
F(i+1,j,2)=rho*u*v;
end 
