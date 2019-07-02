function exp4()

K=3;
M=2;
N=3;
c=0:5:35;
p=10.^(c/10);
test_num=1;
weight=ones(1,K);
eta=1;
noise=0.01;
[~,n]=size(p);
res=zeros(1,n);

for ii=1:1:test_num
    [H_d,H_r,G] = generate_channel1(N,M,K);
    for i=1:1:n
        fprintf('SNR: %i dB ; num: %i \n',c(i),ii);
        res(i)=res(i)+JointFP_CVX(N,M,K,p(i),H_d,H_r,G,weight,eta,noise);
    end
end
res=res/test_num;
plot(0:5:35,res);


end