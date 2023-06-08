%% E1E2 entropy

rng(0);
load('E1E2_1a.mat', 'weight_seq')
load('E1E2_1a.mat', 'msa_aa')

N = size(msa_aa,1);
I = randperm(N);


weight_seq = weight_seq(I);
msa_aa = msa_aa(I,:);

entropy_true = [];
[Profile, Symbols] = seqprofile(msa_aa);
entropy_ind = -Profile.*log(Profile);
S=[];
for i =500:500:N
    [msa_unique ind1 ind2]= unique(msa_aa(1:i,:),'rows');
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
plot(1./[500:500:N],S)

entropy_ind = sum(entropy_ind(~isnan(entropy_ind)));
for chunk = [100:50:1900]
S = [];

for i =chunk:chunk:N
    [msa_unique ind1 ind2]= unique(msa_aa(1:i,:),'rows');
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





p = polyfit(1./[chunk:chunk:N]',S,2);

entropy_true = [entropy_true ;p(3)];
end

%% load inferred 1a E2 model and MCMC samples
load('results\Model_1a.mat','J_MPF_BML')
load('results\MCMC\MCMC_samples_1a_9990.mat', 'samples')

%%  load inferred 1b E2 model and MCMC samples
load('results\Model_1b.mat','J_MPF_BML')
load('results\MCMC\MCMC_samples_1b_9990.mat', 'samples')

%%  load inferred 3a E2 model and MCMC samples
load('results/workspace/test_1226.mat','J_MPF_BML')
load('results\MCMC\MCMC_samples_3a_9990.mat', 'samples')
%% entropy计算

[samples_unique, ind1, ind2]= unique(samples,'rows');


weight_seq = ones(1,size(samples,1))/size(samples,1);
weight_seq_unique=[];
for indi_bin = 1:length(ind1)
    num_term = ind2(ind1(indi_bin));
    ind_values = find(ind2==num_term);
    weight_seq_unique(indi_bin) = sum(weight_seq(ind_values)); %为了得到每个unique sequence占的weight吧
end

Z=0;
E=0;
for i =1:size(samples_unique,1)
    Z = Z+exp(-samples_unique(i,:)*J_MPF_BML*samples_unique(i,:)');
    
    E = E+samples_unique(i,:)*J_MPF_BML*samples_unique(i,:)'*weight_seq_unique(i);
end

entropy_maxent = E+log(Z);
entropy_maxent
%%
E1E2=(entropy_ind - entropy_maxent)/(entropy_ind - max(entropy_true));
I_E1E2 = entropy_ind - max(entropy_true);
E1E2_ind = entropy_ind;
disp(['E1E2']);
disp(['S_ind: ' ,num2str(entropy_ind ),' S_true: ',num2str( max(entropy_true)),' S_maxent: ',num2str( entropy_maxent),' Frac: ',num2str( E1E2)]);


%% E2 entropy

rng(0);
load('E1E2_1a.mat', 'weight_seq')
load('E1E2_1a.mat', 'msa_aa')
range = 193:555;

msa_aa = msa_aa(:,range);

% load('data_fitnessCosts_E2_99900.mat', 'msa_aa')
% load('data_fitnessCosts_E2_99900.mat', 'weight_id')
% weight_seq = weight_id;




N = size(msa_aa,1);
I = randperm(N);

% weight_seq = ones(1,N);

weight_seq = weight_seq(I);
msa_aa = msa_aa(I,:);

entropy_true = [];
[Profile, Symbols] = seqprofile(msa_aa);
entropy_ind = -Profile.*log(Profile);

entropy_ind = sum(entropy_ind(~isnan(entropy_ind)));
for chunk = [100:50:1900]
S = [];

for i =chunk:chunk:N
    [msa_unique ind1 ind2]= unique(msa_aa(1:i,:),'rows');
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





p = polyfit(1./[chunk:chunk:N]',S,2);

entropy_true = [entropy_true ;p(3)];
end

max(entropy_true) ;
load('data_fitnessCosts_E2_99900.mat', 'samples_MCMC')
load('data_fitnessCosts_E2_99900.mat', 'J_MPF_BML')

samples = samples_MCMC;
[samples_unique ind1 ind2]= unique(samples,'rows');


weight_seq = ones(1,size(samples,1))/size(samples,1);
weight_seq_unique=[];
for indi_bin = 1:length(ind1)
    num_term = ind2(ind1(indi_bin));
    ind_values = find(ind2==num_term);
    weight_seq_unique(indi_bin) = sum(weight_seq(ind_values));
end

Z=0;
E=0;
for i =1:size(samples_unique,1)
    Z = Z+exp(-samples_unique(i,:)*J_MPF_BML*samples_unique(i,:)');
    
    E = E+samples_unique(i,:)*J_MPF_BML*samples_unique(i,:)'*weight_seq_unique(i);
end

entropy_maxent = E+log(Z);
E2=(entropy_ind - entropy_maxent)/(entropy_ind - max(entropy_true)) ;

E2_max = entropy_maxent;
disp(['E2']);
disp(['S_ind: ' ,num2str(entropy_ind ),' S_true: ',num2str( max(entropy_true)),' S_maxent: ',num2str( entropy_maxent),' Frac: ',num2str( E2)]);
%% E2_new entropy

rng(0);
load('E2_1a.mat', 'weight_seq')
load('E2_1a.mat', 'msa_aa')
% range = 193:555;
% 
% msa_aa = msa_aa(:,range);

% load('data_fitnessCosts_E2_99900.mat', 'msa_aa')
% load('data_fitnessCosts_E2_99900.mat', 'weight_id')
% weight_seq = weight_id;




N = size(msa_aa,1);
I = randperm(N);

% weight_seq = ones(1,N);

weight_seq = weight_seq(I);
msa_aa = msa_aa(I,:);

entropy_true = [];
[Profile, Symbols] = seqprofile(msa_aa);
entropy_ind = -Profile.*log(Profile);

entropy_ind = sum(entropy_ind(~isnan(entropy_ind)));
for chunk = [100:50:1900]
S = [];

for i =chunk:chunk:N
    [msa_unique ind1 ind2]= unique(msa_aa(1:i,:),'rows');
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





p = polyfit(1./[chunk:chunk:N]',S,2);

entropy_true = [entropy_true ;p(3)];
end

max(entropy_true) ;
load('E2_1a.mat', 'samples')
load('E2_1a.mat', 'J_MPF_BML')

% samples = samples_MCMC;
[samples_unique ind1 ind2]= unique(samples,'rows');


weight_seq = ones(1,size(samples,1))/size(samples,1);
weight_seq_unique=[];
for indi_bin = 1:length(ind1)
    num_term = ind2(ind1(indi_bin));
    ind_values = find(ind2==num_term);
    weight_seq_unique(indi_bin) = sum(weight_seq(ind_values));
end

Z=0;
E=0;
for i =1:size(samples_unique,1)
    Z = Z+exp(-samples_unique(i,:)*J_MPF_BML*samples_unique(i,:)');
    
    E = E+samples_unique(i,:)*J_MPF_BML*samples_unique(i,:)'*weight_seq_unique(i);
end

entropy_maxent = E+log(Z);
E2_new=(entropy_ind - entropy_maxent)/(entropy_ind - max(entropy_true)) ;

E2_max_new = entropy_maxent;
disp(['E2_new']);
disp(['S_ind: ' ,num2str(entropy_ind ),' S_true: ',num2str( max(entropy_true)),' S_maxent: ',num2str( entropy_maxent),' Frac: ',num2str( E2_new)]);
%% E1 entropy

rng(0);
load('E1E2_1a.mat', 'weight_seq')
load('E1E2_1a.mat', 'msa_aa')
range = 1:192;

msa_aa = msa_aa(:,range);
load('samples_E1.mat')
samples = all_samples{1};




N = size(msa_aa,1);
I = randperm(N);

% weight_seq = ones(1,N);

weight_seq = weight_seq(I);
msa_aa = msa_aa(I,:);

entropy_true = [];
[Profile, Symbols] = seqprofile(msa_aa);
entropy_ind = -Profile.*log(Profile);

entropy_ind = sum(entropy_ind(~isnan(entropy_ind)));
for chunk = [100:50:1900]
S = [];

for i =chunk:chunk:N
    [msa_unique ind1 ind2]= unique(msa_aa(1:i,:),'rows');
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





p = polyfit(1./[chunk:chunk:N]',S,2);

entropy_true = [entropy_true ;p(3)];
end


load('E1_1a.mat', 'J_MPF_BML')

[samples_unique ind1 ind2]= unique(samples,'rows');


weight_seq = ones(1,size(samples,1))/size(samples,1);
weight_seq_unique=[];
for indi_bin = 1:length(ind1)
    num_term = ind2(ind1(indi_bin));
    ind_values = find(ind2==num_term);
    weight_seq_unique(indi_bin) = sum(weight_seq(ind_values));
end

Z=0;
E=0;
for i =1:size(samples_unique,1)
    Z = Z+exp(-samples_unique(i,:)*J_MPF_BML*samples_unique(i,:)');
    
    E = E+samples_unique(i,:)*J_MPF_BML*samples_unique(i,:)'*weight_seq_unique(i);
end

entropy_maxent = E+log(Z);

E1=(entropy_ind - entropy_maxent)/(entropy_ind - max(entropy_true)) ;
E1_max = entropy_maxent;
disp(['E1']);
disp(['S_ind: ' ,num2str(entropy_ind ),' S_true: ',num2str( max(entropy_true)),' S_maxent: ',num2str( entropy_maxent),' Frac: ',num2str( E1)]);


disp(['Independent']);

indept =(E1E2_ind-E1_max-E2_max)/I_E1E2;
disp([' S_maxent: ',num2str( E1_max+E2_max),' Frac: ',num2str( indept)]);


disp(['Independent new']);

indept =(E1E2_ind-E1_max-E2_max_new)/I_E1E2;
disp([' S_maxent: ',num2str( E1_max+E2_max_new),' Frac: ',num2str( indept)]);
%%
run startup.m

% old 67.4
Entropy = [103.2 40.9 62.9 7.9];
% Entropy = [103.2 40.9 67.4 7.9];
FIG=figure;


hold on;

% l = categorical({'1a','1b'});
% l = reordercats(l,{'1a','1b'});
% b = bar(l,[length(cum_fre_1a) length(cum_fre_1b)],0.2,'LineWidth',0.1,'BarWidth',1);
% b.FaceColor = 'flat';
% b.CData(1,:)=purple;
bar(1,Entropy(1), 0.8, 'FaceColor',color_scheme(9,:),'LineWidth',0.2,'FaceAlpha',0.6,'EdgeColor','none');

bar(3,Entropy(2), 0.8, 'FaceColor',color_scheme(1,:),'LineWidth',0.2,'FaceAlpha',0.6,'EdgeColor','none');

bar(2,Entropy(3), 0.8, 'FaceColor',color_scheme(2,:),'LineWidth',0.2,'FaceAlpha',0.6,'EdgeColor','none');

bar(4,Entropy(4), 0.8, 'FaceColor',color_scheme(4,:),'LineWidth',0.2,'FaceAlpha',0.6,'EdgeColor','none');

% b.CData(2,:)=orange;
% for i =1:4
%     bar(i,Entropy(i), 0.8, 'FaceColor',color_scheme(i,:),'LineWidth',0.2);
% end
xlim([0.5 4.5])
set(gca,'XTick',1:4,'XTickLabel',...
    {'S_{ind}','S_{joint}','S_{individual}','S_{true}'});
set(gca,'TickDir','out')
FIG.Name = 'peak_tract'
set(gca,'TickLength',[0.02, 0.01])
ylabel({'Entropy (nats)'})
% legend('Subtype 1a', 'Subtype 1b','Location','best')
FIG.Units = 'centimeters';


set(gcf,'Position',[10 10 8 10]);
set(gca,'Position',[.18 .17 .8 .74]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'Position',[.22 .25 .88 .71]);  %调整 XLABLE和YLABLE不会被切掉
% set(gcf,'Position',[10 10 7.84 6]);
% set(gca,'Position',[.18 .17 .76 .74]);  %调整 XLABLE和YLABLE不会被切掉

figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
ylim([0 110])
set(gca,'YTick', [0 50 100])
set(gca,'TickLength',[0.035, 0.03])
try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end
%%
% new
Frac = [ 0.65407 0.42321];


% old
% Frac = [ 0.65407 0.37611];


FIG=figure;


hold on;

% l = categorical({'1a','1b'});
% l = reordercats(l,{'1a','1b'});
% b = bar(l,[length(cum_fre_1a) length(cum_fre_1b)],0.2,'LineWidth',0.1,'BarWidth',1);
% b.FaceColor = 'flat';
% b.CData(1,:)=purple;
% b.CData(2,:)=orange;

bar(1,Frac(1), 0.5, 'FaceColor',color_scheme(1,:),'LineWidth',0.2,'FaceAlpha',0.6,'EdgeColor','none');

bar(2,Frac(2), 0.5, 'FaceColor',color_scheme(2,:),'LineWidth',0.2,'FaceAlpha',0.6,'EdgeColor','none');
xlim([0.5 2.5])
set(gca,'XTick',1:2,'XTickLabel',...
    {'JM','IM'});
set(gca,'TickDir','out')
FIG.Name = 'peak_tract'
set(gca,'TickLength',[0.02, 0.01])
% ylabel({'I_{model}/I'})
ylabel({'Fraction of the correlated structure (FCS)', 'captured by the model'})

set(gca,'YTick',0:0.2:1)

% legend('Subtype 1a', 'Subtype 1b','Location','best')
FIG.Units = 'centimeters';


set(gcf,'Position',[10 10 6 8]);
set(gca,'Position',[.24 .17 .72 .65]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'Position',[.22 .25 .88 .71]);  %调整 XLABLE和YLABLE不会被切掉
% set(gcf,'Position',[10 10 7.84 6]);
% set(gca,'Position',[.18 .17 .76 .74]);  %调整 XLABLE和YLABLE不会被切掉

figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.24, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
ylim([0 1])
set(gca,'YTick', [0 0.5 1])
set(gca,'TickLength',[0.035, 0.03])
try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end