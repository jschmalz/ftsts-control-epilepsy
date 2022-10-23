function [W_EE, W_EI, W_IE, W_II] = make_weights(sigma_E,sigma_I,N_E)

% weight matrix
sigma = ceil(sigma_E*N_E);
W_EE = zeros(N_E,N_E);

for i = 1:N_E
   for j = 1:N_E
      W_EE(i,j) = mynormpdf(abs(i-j),0,sigma); 
   end
end
W_EE = (104)*W_EE./sum(W_EE,1);

sigma = ceil(sigma_E*N_E);
W_EI = zeros(N_E,N_E);



for i = 1:N_E
   for j = 1:N_E
      W_EI(i,j) = mynormpdf(abs(i-j),0,sigma); 
   end
end
W_EI = (100)*W_EI./sum(W_EI,1);

sigma = ceil(sigma_I*N_E);
W_IE = zeros(N_E,N_E);

for i = 1:N_E
   for j = 1:N_E
      W_IE(i,j) = mynormpdf(abs(i-j),0,sigma); 
   end
end

W_IE = 250*W_IE./sum(W_IE,1) + 50*(1/N_E);

sigma = ceil(sigma_I*N_E);
W_II = zeros(N_E,N_E);

for i = 1:N_E
   for j = 1:N_E
      W_II(i,j) = mynormpdf(abs(i-j),0,sigma); 
   end
end

W_II =  250*W_II./sum(W_II,1) + 50*(1/N_E);

end