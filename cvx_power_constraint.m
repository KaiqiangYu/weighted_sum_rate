function res=cvx_power_constraint(lambda_x,H_k,beta,constant_alpha,M,K)

temp_a=zeros(M,M);
for i=1:1:K
    temp_a=temp_a+abs(beta(i))^2*H_k(:,i)*H_k(:,i)';
end
temp_c=diag(lambda_x*ones(1,M));
temp_b=(temp_a+temp_c)^-1;
W=zeros(M,K);

for i=1:1:K
      W(:,i)=sqrt(constant_alpha(i))*beta(i)*temp_b*H_k(:,i);
end

res=0;
for i=1:1:K
     res=res+norm(W(:,i))^2;
end

end