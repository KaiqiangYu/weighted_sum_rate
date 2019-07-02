function [H_d,H_r,G] = generate_channel(N,M,K)
%Generate system model consisting of K communication pairs and L IRS
%   K: the number of commnunication pairs(source-destination pairs with signal antenna)
%   L: the number of intelligence reflecting surface (IRS)
%   M_l: the number of elements of each IRS
%   @author Hejl Jun,2019&Supported by Yui
%   @copyright Shanghaitech&KIAA IC1011 YZLHX
%
%%
T=1000;
alpha=2.3;
circle_radius=50;
distance=300;
source_position=[0,0];
destination_position=zeros(K,2);
for i=1:1:K
    y=2*pi*rand;
    destination_position(i,1)=distance+circle_radius*cos(y);
    destination_position(i,2)=circle_radius*sin(y);
end

%IRS position
IRS=[100,50];

S_D=zeros(1,K);
I_D=zeros(1,K);
S_I=sqrt(T*norm(source_position-IRS)^-alpha);
for i=1:1:K
    S_D(i)=sqrt(T*norm(destination_position(i,:))^-alpha);
    I_D(i)=sqrt(T*norm(destination_position(i,:)-IRS)^-alpha);
end


H_d=zeros(M,K);
H_r=zeros(N,K);
for i=1:1:K
    H_d(:,i)=S_D(i)*(normrnd(0,1/sqrt(2),M,1)+1i* normrnd(0,1/sqrt(2),M,1));
    H_r(:,i)=I_D(i)*(normrnd(0,1/sqrt(2),N,1)+1i* normrnd(0,1/sqrt(2),N,1));
end


G=S_I*(normrnd(0,1/sqrt(2),N,M)+1i* normrnd(0,1/sqrt(2),N,M));



end