%% estimate the true entropy
load('Model_1a.mat', 'weight_id')
load('Model_1a.mat', 'msa_aa')

weight_seq = weight_id;
N = size(msa_aa,1);

all_S = [];
num_iter = 100;

%for i =[500:500:N N]
for i =[N]
    S=[];
    for j =1:num_iter
        rng(j+i);
        N = size(msa_aa,1);
        I = randperm(N);
        
        
        weight_seq = weight_seq(I);
        msa_aa_tmp = msa_aa(I,:);
        [msa_unique ind1 ind2]= unique(msa_aa_tmp(1:i,:),'rows');
        weight_seq_unique=[];
        for indi_bin = 1:length(ind1)
            num_term = ind2(ind1(indi_bin));
            ind_values = find(ind2==num_term);
            weight_seq_unique(indi_bin) = sum(weight_seq(ind_values));
        end
        weight_seq_unique = weight_seq_unique/sum(weight_seq_unique);
        
        Stmp = sum(-weight_seq_unique.*log(weight_seq_unique));
        S = [S;Stmp];

    end
    all_S =[all_S S];
end
%%
run function/startup.m
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

p = polyfit(1./[500:500:N N],mean(all_S,1),2);
FIG = figure;
scatter(1./[500:500:N N]',mean(all_S,1), 20,blue,'filled','MarkerFaceAlpha',0.6); 
xlabel('1/M')
ylabel('\langle S_{naive} \rangle')

x = 0:max(1./[500:500:N N])/100:max(1./[500:500:N N]);
y = x.^2*p(1)+x*p(2)+p(3);

hold on ;

plot(x,y,'--','Color',blue)

legend('\langle S_{naive} \rangle','fit line')
FIG.Units = 'centimeters';
set(gca,'Position',[.15 .2 .8 .77]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'Position',[.1 .13 .2 .83]);  %调整 XLABLE和YLABLE不会被切掉

set(gcf,'Position',[10 10 8 6]);
% set(gca,'Position',[.1 .13 .67 .83]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(gca,'TickDir','out')
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(gca,'Position',[.1 .12 .74 .87]);  %调整 XLABLE和YLABLE不会被切掉
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.145, 0.5, 0]);
set(get(gca,'XLabel'), 'Units', 'Normalized', 'Position', [0.5, -0.16, 0]);
set(get(gca,'title'), 'Units', 'Normalized', 'Position', [0.42, 1.1, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'TickLength',[0.02, 0.03])

%% True entropy 1a
load('Model_1a.mat', 'weight_id')
load('Model_1a.mat', 'msa_aa')

weight_seq = weight_id;
num_seqs = size(msa_aa,1);
L=363;
M_list=[];
S_list=[];

% subsampling
for M=[1:50:1121 1121]
    S_M=[];
    for repeat = 1:5
        rng(M+repeat);
        S=0;
        index = randperm(num_seqs,M)';
        population = msa_aa(index,:);
        population_weight = weight_seq(index);
        unique_seq = unique(population,'rows');
        N= sum(population_weight);

        for i = 1:size(unique_seq,1)
            seqi = unique_seq(i,:);
            flag = sum(population==seqi,2)==L;
            freq_seqi = sum(population_weight(flag))/N;
            S = S - freq_seqi*log(freq_seqi);
        end
        S_M = [S_M;S];
        
    end
    M_list=[M_list;M];
    S_list=[S_list;mean(S_M)];
end

x_1a= M_list;
y_1a= S_list;
scatter(x_1a,y_1a);

p_1a = polyfit(x_1a,y_1a,3); %% 目前觉得3比较吻合
x1=500:100:10000;
y1= polyval(p_1a,x1);

hold on
plot(x1,y1)

xlabel('Number of sequences')
ylabel('Ture Entropy')
legend("True entropy","estimated curve for true entropy")
title("Fitting curve for the true entropy 【1a E2 protein】")
%% 
load('Model_1b.mat', 'weight_seq')
load('Model_1b.mat', 'msa_aa')

num_seqs = size(msa_aa,1);
L=320; % 1b这里的msa是去掉了conserved aa的，不过也不影响我们算weighted seq probability
M_list=[];
S_list=[];

% subsampling
for M=[1:50:1121 1121]
    S_M=[];
    for repeat = 1:5
        rng(M+repeat);
        S=0;
        index = randperm(num_seqs,M)';
        population = msa_aa(index,:);
        population_weight = weight_seq(index);
        unique_seq = unique(population,'rows');
        N= sum(population_weight);

        for i = 1:size(unique_seq,1)
            seqi = unique_seq(i,:);
            flag = sum(population==seqi,2)==L;
            freq_seqi = sum(population_weight(flag))/N;
            S = S - freq_seqi*log(freq_seqi);
        end
        S_M = [S_M;S];
        
    end
    M_list=[M_list;M];
    S_list=[S_list;mean(S_M)];
end

x_1b= M_list;
y_1b= S_list;
scatter(x_1b,y_1b);

p_1b = polyfit(x_1b,y_1b,3); %% 目前觉得3比较吻合
x2=500:100:10000;
y2= polyval(p_1b,x2);

hold on
plot(x2,y2)

xlabel('Number of sequences')
ylabel('Ture Entropy')
legend("True entropy","estimated curve for true entropy")
title("Fitting curve for the true entropy 【1b E2 protein】")
%% 3a E2
load('results/workspace/test_1226.mat', 'weight_seq')
load('results/workspace/test_1226.mat', 'msa_aa')

num_seqs = size(msa_aa,1);
L=340; % 1b这里的msa是去掉了conserved aa的，不过也不影响我们算weighted seq probability
M_list=[];
S_list=[];

% subsampling
for M=[1:50:1121 1121]
    S_M=[];
    for repeat = 1:5
        rng(M+repeat);
        S=0;
        index = randperm(num_seqs,M)';
        population = msa_aa(index,:);
        population_weight = weight_seq(index);
        unique_seq = unique(population,'rows');
        N= sum(population_weight);

        for i = 1:size(unique_seq,1)
            seqi = unique_seq(i,:);
            flag = sum(population==seqi,2)==L;
            freq_seqi = sum(population_weight(flag))/N;
            S = S - freq_seqi*log(freq_seqi);
        end
        S_M = [S_M;S];
        
    end
    M_list=[M_list;M];
    S_list=[S_list;mean(S_M)];
end

x_3a= M_list;
y_3a= S_list;
scatter(x_3a,y_3a);

p_3a = polyfit(x_3a,y_3a,3); %% 目前觉得3比较吻合
x3=500:100:10000;
y3= polyval(p_3a,x3);

hold on
plot(x3,y3)

xlabel('Number of sequences')
ylabel('Ture Entropy')
legend("True entropy","estimated curve for true entropy")
title("Fitting curve for the true entropy 【3a E2 protein】")

%% plot 3 fitting curves
% 1a
a=-4.2965e+03; b=-1.9565e-04; c=4.2965e+03;
%x = 1000:100:1e10;
%y = a.* x .^b+c;
%plot(x,y)

syms M
f = a*M^b + c;
S_true_1a = limit(f,M,Inf)
%% 1b
a=44.4298; b=0.0195; c=-44.4298; %这个b的是正数哎！！
%x = 1000:100:1e10;
%y = a.* x .^b+c;
%plot(x,y)

syms M
f = a*M^b + c;
S_true_1a = limit(f,M,Inf)

%% 3a


%% MCMC model entropy
load("results/MCMC/MCMC1e7_samples_1a.mat")
load("results/Model_1a.mat","J_MPF_BML")

E = diag(samples*J_MPF_BML*samples'); %all
[unique_seq, ind1]=unique(samples,'rows');
unique_seq_weight=zeros(size(unique_seq,1),1);

for i = 1:size(unique_seq,1)
    seqi= unique_seq(i,:);
    num = sum(sum(samples==seqi,2)==624);
    unique_seq_weight(i) = num/size(unique_seq,1);
end
%% 

% case 1: 0.9376
Q = exp(-E);
Q = Q ./ sum(Q);
unique_Q = Q(ind1);
freq = unique_Q .* unique_seq_weight;
prob_seq = freq ./ sum(freq);
S = sum(-prob_seq .* log(prob_seq));
%% 

% case 2: 1.1664
Q = exp(-E);
Q = Q ./ sum(Q);
unique_Q = Q(ind1);
freq = unique_Q;
prob_seq = freq ./ sum(freq);
S = sum(-prob_seq .* log(prob_seq));

%% case 3: 0.9376
Q = exp(-E(ind1));
Q = Q ./ sum(Q);
freq = Q .* unique_seq_weight;
prob_seq = freq ./ sum(freq);
S = sum(-prob_seq .* log(prob_seq));

%% case 4:1.1664
Q = exp(-E(ind1));
freq = Q ./ sum(Q)
prob_seq = freq;
S = sum(-prob_seq .* log(prob_seq));

%% 按照locationa看Model entropy
load("results/MCMC/MCMC1e7_samples_1a.mat")
load("results/Model_1a.mat","J_MPF_BML")















