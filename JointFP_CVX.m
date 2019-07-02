function res=JointFP_CVX(N,M,K,p_max,H_d,H_r,G,weight,eta,noise)
%% Fractional programming for wighted Sum-Rate Optimization Problem
%   K: the number of single-antenna users
%   N: the number of antennas of AP
%   M: the number of passive elements in IRS
%   @author Hejl Jun,2019&Supported by Yui
%   @copyright Shanghaitech&KIAA IC1011 YZLHX
%   @thanks 5-307, Software College, SDU & 1C-401, SIST, SHT
%
%% initialize the variables at the beginning of FP procedures
% theta: reflection coefficient vector
% W: beamforming matrix
theta=ones(N,1);
theta=theta*exp(-j);

W = normrnd(0,1/sqrt(2),M,K)+1i* normrnd(0,1/sqrt(2),M,K);

%% aggregate different channels as one vector form
% h_k=h_{d,k}+G^H*theta*h_k



%initialize the variable
res_last=30;
res=0;
count=0;

%% FP procedures to solve multi-ratio fractional programming

while abs(res-res_last)>0.0001

H_k=zeros(M,K);
for i=1:1:K
    H_k(:,i)=H_d(:,i)+G'*diag(theta)*H_r(:,i);
end
    
count=count+1;
res_last=res;
% for fixed theta and W, compute alpha
alpha=zeros(1,K);
for i=1:1:K
      temp_a=abs(H_k(:,i)'*W(:,i))^2;  %compute the numerator
      temp_b=0;
      for m=1:1:K
           if m~=i
                 temp_b=temp_b+abs(H_k(:,m)'*W(:,m))^2; %compute the denominator
           end
      end
      temp_b=temp_b+noise;
      alpha(i)=(temp_a)/(temp_b); %combine the above
end

% for the fixed theta and W, and updated alpha, compute constant_alpha
constant_alpha=zeros(1,K);
for i=1:1:K
      constant_alpha(i)=weight(i)*(1+alpha(i));
end

% for the fixed theta and constant_alpha, optimize beamforming vector W
% firstly, for fixed W, compute beta
beta=zeros(1,K);
for i=1:1:K
      temp_a=sqrt(constant_alpha(i))*(H_k(:,i)'*W(:,i));  %compute the numerator
      temp_b=0;
      for m=1:1:K
            temp_b=temp_b+abs(H_k(:,m)'*W(:,m))^2; %compute the denominator
      end
      temp_b=temp_b+noise;
      beta(i)=(temp_a)/(temp_b); %combine the above
end

% for the fixed beta, compute W
% firstly compute the optimal dual variable lambda
% cvx_begin
% variable lambda_x(1,1)
% expression aaaa
% aaa=1/(lambda_x*diag(ones(1,3)))
% minimize lambda_x
% subject to
%    lambda_x>=0
%    cvx_power_constraint(lambda_x,H_k,beta,constant_alpha,M,K)<=p_max
% cvx_end
lambda_up=1000;
lambda_down=0;
lambda_x=0;
while abs(lambda_up-lambda_down)>0.0001
    lambda_x=(lambda_up+lambda_down)/2;
    if cvx_power_constraint(lambda_x,H_k,beta,constant_alpha,M,K)>p_max
             lambda_down=lambda_x;
    else
             lambda_up=lambda_x;
    end
end

% compute W
temp_wa=zeros(M,M);
for i=1:1:K
    temp_wa=temp_wa+abs(beta(i))^2*H_k(:,i)*H_k(:,i)';
end
temp_wb=(lambda_x*diag(ones(1,M))+temp_wa)^-1;
for i=1:1:K
      W(:,i)=sqrt(constant_alpha(i))*beta(i)*temp_wb*H_k(:,i);
end

%  for fixed W, optimize theta via two different methods (CVX and DC)
%  inital other variables
for i=1:1:K
     temp_var_a=zeros(N,K);
     var_b=zeros(K,K);
     for m=1:1:K
          temp_var_a(:,m)=sqrt(eta)*diag(H_r(:,m)')*G*W(:,i);
          var_b(i,m)=H_d(:,m)'*W(:,i);
     end
     var_a{i}=temp_var_a;
end

%  for fixed W and theta, compute eps
eps=zeros(1,K);
for i=1:1:K
      temp_a=sqrt(constant_alpha(i))*(var_b(i,i)+theta'*var_a{i}(:,i));  %compute the numerator
      temp_b=0;
      for m=1:1:K
            temp_b=temp_b+abs(var_b(m,i)+theta'*var_a{m}(:,i))^2; %compute the denominator
      end
      temp_b=temp_b+noise;
      eps(i)=(temp_a)/(temp_b); %combine the above
end

% for fixed eps, compute theta by CVX
cvx_begin quiet
variable var_theta(N,1) complex
maximize obj_fun(constant_alpha,eps,var_theta,var_a,var_b,noise,K)
subject to
        var_theta.*conj(var_theta)<=ones(N,1)
cvx_end

theta=var_theta;
% projection to non-convex constraint

% for i=1:1:N
%     cvx_begin quiet
%     variable vv_theta(1,1) complex 
%     minimize (vv_theta-theta(i))'*(vv_theta-theta(i))
%     subject to
%             vv_theta*vv_theta'==1
%     cvx_end
%     
%     theta(i)=vv_theta;
% end

temp_res=0;
for i=1:1:K
     temp_a=abs(H_k(:,i)'*W(:,i))^2;  %compute the numerator
      temp_b=0;
      for m=1:1:K
           if m~=i
                 temp_b=temp_b+abs(H_k(:,m)'*W(:,m))^2; %compute the denominator
           end
      end
      temp_b=temp_b+noise;
      temp_res=temp_res+weight(i)*log(1+(temp_a)/(temp_b))/log(2); %combine the above
    
end
res=temp_res;

fprintf('obj: %d  bound: %d \n',res,abs(res-res_last))

if count>=300
    fprintf('not converence!!! error: %d related error: %d',abs(res-res_last),abs(res-res_last)/res);
    break;
end


end

    



end
