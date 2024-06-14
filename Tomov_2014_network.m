 clear all
H = 0
g_ex = 0.6; g_in = 1;

RS_p  = 0.64; RS_val  = [0.02 0.2  -65 8];
IB_p  = 0   ; IB_val  = [0.02 0.2  -55 4];
CH_p  = 0.16; CH_val  = [0.02 0.2  -50 2];
FS_p  = 0  ;  FS_val  = [0.1  0.2  -65 2];
LTS_p = 0.2;  LTS_val = [0.02 0.25 -65 2];

if H == 0; N = 2^9;
elseif H == 1; N = 2^9; %%
end
N_ex = round((RS_p+IB_p+CH_p)*N);
N_in = round((FS_p+LTS_p)*N);

A = [repmat(RS_val(1),round(RS_p*N),1); repmat(IB_val(1),round(IB_p*N),1); repmat(CH_val(1),round(CH_p*N),1) ...
    repmat(FS_val(1),round(FS_p*N),1); repmat(LTS_val(1),round(LTS_p*N),1)] ;
B = [repmat(RS_val(2),round(RS_p*N),1); repmat(IB_val(2),round(IB_p*N),1); repmat(CH_val(2),round(CH_p*N),1) ...
    repmat(FS_val(2),round(FS_p*N),1); repmat(LTS_val(2),round(LTS_p*N),1)] ;
C = [repmat(RS_val(3),round(RS_p*N),1); repmat(IB_val(3),round(IB_p*N),1); repmat(CH_val(3),round(CH_p*N),1) ...
    repmat(FS_val(3),round(FS_p*N),1); repmat(LTS_val(3),round(LTS_p*N),1)] ;
D = [repmat(RS_val(4),round(RS_p*N),1); repmat(IB_val(4),round(IB_p*N),1); repmat(CH_val(4),round(CH_p*N),1) ...
    repmat(FS_val(4),round(FS_p*N),1); repmat(LTS_val(4),round(LTS_p*N),1)] ;
tmp = randperm(N_ex); A(1:N_ex) = A(tmp,1); B(1:N_ex) = B(tmp,1); C(1:N_ex) = C(tmp,1); D(1:N_ex) = D(tmp,1);
tmp = randperm(N_in)+N_ex; A(N_ex+1:N) = A(tmp,1); B(N_ex+1:N) = B(tmp,1); C(N_ex+1:N) = C(tmp,1); D(N_ex+1:N) = D(tmp,1);

if H == 0; p = 0.02;
    connect = zeros(N,N);
    connect(randperm(N*N,round(N*N*p))) = 1;
    
    g_coef = [repmat(g_ex,N_ex,1); repmat(g_in,N_in,1)];
    G = [connect(:,1:N_ex)*g_ex connect(:,(N_ex+1):N)*g_in]; %ini
    E = [repmat(0,N_ex,1); repmat(-80,N_in,1)];
    tau = [repmat(5,N_ex,1); repmat(6,N_in,1)];
elseif H == 1; p = 0.01; %%
    connect_ex1_mod = zeros(N/2,floor(N_ex/2));
    connect_ex1_int = zeros(N/2,floor(N_ex/2));
    tmp = numel(connect_ex1_mod);
    connect_ex1_mod(randperm(tmp,round(tmp*p*(2-0.1)))) = 1;
    connect_ex1_int(randperm(tmp,round(tmp*p*(0.1)))) = 1;
    
    connect_ex2_int = zeros(N/2,ceil(N_ex/2));
    connect_ex2_mod = zeros(N/2,ceil(N_ex/2));
    tmp = numel(connect_ex2_mod);
    connect_ex2_mod(randperm(tmp,round(tmp*p*(2-0.1)))) = 1;
    connect_ex2_int(randperm(tmp,round(tmp*p*(0.1)))) = 1;
    
    connect_in1_mod = zeros(N/2,ceil(N_in/2));
    connect_in1_int = zeros(N/2,ceil(N_in/2));
    tmp = numel(connect_in1_mod);
    connect_in1_mod(randperm(tmp,round(tmp*p*(2)))) = 1;
    connect_in1_int(randperm(tmp,round(tmp*p*(0)))) = 1;
    
    connect_in2_int = zeros(N/2,floor(N_in/2));
    connect_in2_mod = zeros(N/2,floor(N_in/2));
    tmp = numel(connect_in2_mod);
    connect_in2_mod(randperm(tmp,round(tmp*p*(2)))) = 1;
    connect_in2_int(randperm(tmp,round(tmp*p*(0)))) = 1;
    
    connect = [[connect_ex1_mod; connect_ex1_int] [connect_ex2_int; connect_ex2_mod] ...
        [connect_in1_mod; connect_in1_int] [connect_in2_int; connect_in2_mod]];
   
    g_coef = [repmat(g_ex,N_ex,1); repmat(g_in,N_in,1)];
    G = [connect(:,1:N_ex)*g_ex connect(:,(N_ex+1):N)*g_in]; %ini
    E = [repmat(0,N_ex,1); repmat(-80,N_in,1)];
    tau = [repmat(5,N_ex,1); repmat(6,N_in,1)];
    
    n_ex1 = size(connect_ex1_mod,2);
    n_ex2 = size(connect_ex2_mod,2);
    n_in1 = size(connect_in1_mod,2);
    n_in2 = size(connect_in2_mod,2);
    H1_sort = [1:n_ex1 N_ex+1:N_ex+n_in1 n_ex1+1:N_ex N-n_in2+1:N];
    connect = connect(:,H1_sort);
    
    A = A(H1_sort); B = B(H1_sort); C = C(H1_sort); D = D(H1_sort);
    G = G(:,H1_sort); g_coef = g_coef(H1_sort); E = E(H1_sort); tau = tau(H1_sort); 
end
% figure
% subplot(2,2,1)
% imagesc(connect)
% subplot(2,2,2)
% imagesc(B)
% subplot(2,2,3)
% imagesc(G)
% subplot(2,2,4)
% imagesc(tau)


h = 0.01;
t = 1:h:100; firings = []; 
v_tmp = NaN(N,length(t)); u_tmp = NaN(N,length(t)); 
G_tmp = NaN(N,length(t)); I_tmp = NaN(N,length(t));
v = C+repmat(80,N,1).*rand(N,1);                         %ini  
u = repmat(0,N,1)+repmat(10,N,1).*rand(N,1);             %ini
v_dot =@(v,u,I)(0.04*(v.^2)+5*v+140-u+I);
u_dot =@(v,u)(A.*(B.*v-u));
G_dot =@(G,tau)(-G./(tau'))
for i = 1:(length(t)-1)
    
    %
    fired=find(v>=30); not_fired=find(~(v>=30));
    firings=[firings; t(i)+0*fired,fired];
    v(fired)=C(fired);
    u(fired)=u(fired)+D(fired);

    G(:,fired) = G(:,fired)+ connect(:,fired).*g_coef(fired)';
    I = zeros(N,1) + sum(G.*(repmat(E',N,1)-connect.*v),2);
    if t(i) < 100
        P_stim = zeros(N,1);
        P_stim(randperm(N,N/2)) = 1;
        %P_stim(randperm(N,N)) = 1;
        I = I + repmat(10,N,1).*P_stim + repmat(10,N,1).*P_stim.*rand(N,1);
    end
    G(:,not_fired) = G(:,not_fired) + h*(-G(:,not_fired)./(tau(not_fired)'));
%     j1 = G_dot(G(:,not_fired),tau(not_fired));
%     j2 = G_dot(G(:,not_fired)+h*connect(:,not_fired),tau(not_fired));
%     G(:,not_fired) = G(:,not_fired)+(0.5*j1+0.5*j2)*h;
    
    %RK2 Heun's method
    k1 = v_dot(v,u,I);
    g1 = u_dot(v,u);
    k2 = v_dot(v+h,u+h,I);
    g2 = u_dot(v+h,u+h);
    v = v+(0.5*k1+0.5*k2)*h;
    u = u+(0.5*g1+0.5*g2)*h;
    
    if mod(i, 3) == 0
        disp([num2str(t(i),'%4.2f') 'ms ' num2str(length(firings)) ' fired'])
    end
    v_tmp(:,i) = v;
    u_tmp(:,i) = u;
    G_tmp(:,i) = sum(G,2);
    I_tmp(:,i) = I;
end

% cd E:\20191019_Michigan_Comp_NSC\Mich_Spike_model
% save(['Tomov_N' num2str(N) '_p' strrep(num2str(p),'.','') '_gex' strrep(num2str(g_ex),'.','') '_gin' strrep(num2str(g_in),'.','')],'t','v_tmp','u_tmp','G_tmp', 'firings')

neuron = [1 50 500]
figure('Renderer', 'painters', 'Position', [10 10 1500 800])
subplot(5,1,1)
plot(firings(:,1),firings(:,2),'.'); xlim([1 max(t)]); 
subplot(5,1,2)
plot(t,v_tmp(neuron,:)); xlim([1 max(t)]); ylabel('v')
subplot(5,1,3)
plot(t,u_tmp(neuron,:)); xlim([1 max(t)]); ylabel('u')
subplot(5,1,4)
plot(t,G_tmp(neuron,:)); xlim([1 max(t)]); ylabel('mean(G_{axon})')
subplot(5,1,5)
plot(t,I_tmp(neuron,:)); xlim([1 max(t)]); ylabel('I')


