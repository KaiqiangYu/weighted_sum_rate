function res=obj_fun(constant_alpha,eps,var_theta,var_a,var_b,noise,K)

temp_a=0;
for i=1:1:K
    temp_a=temp_a+2*sqrt(constant_alpha(i))*real(conj(eps(i))*var_theta'*var_a{i}(:,i)+conj(eps(i))*var_b(i,i));
end

temp_b=0;
for i=1:1:K
    temp_c=0;
    for m=1:1:K
       % (var_b(m,i)+var_theta'*var_a{m}(:,i))*(var_b(m,i)+var_theta'*var_a{m}(:,i))'
          temp_c=temp_c+(var_b(m,i)+var_theta'*var_a{m}(:,i))*(var_b(m,i)+var_theta'*var_a{m}(:,i))';
    end
    temp_b=temp_b+abs(eps(i))^2*(temp_c+noise);
end

res=temp_a-temp_b;



end