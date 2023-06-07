%% load 1a model
load('results\Model_1a.mat')

%% load 1b model
load('results\Model_1b.mat')

%% load inferred 3a E2 model 
load('results/workspace/test_1226.mat')
%% 
rng(0)
rand

%% Metropolis MC sampling
rng(0)
stream = RandStream("mcg16807","Seed",0);
tic;
samples=[];
N=1e7; %1e7
Burn_in=0; %1e4
Thin_in=1; %1e3
USE_H=0;

num_mutants_combine_array_acc = cumsum(num_mutants_combine_array);
num_mutants_combine_array_acc_all = [0 num_mutants_combine_array_acc];
positions = size(num_mutants_combine_array,2);    
msa_bin_unique = unique(msa_bin,'row');
randvalue = randi(stream,[1 size(msa_bin_unique,1)]);
x = msa_bin_unique(randvalue,:);      % randomly pick a sequence from unique MSA sequence
n=0;

if(USE_H)
    energy = x*H';
else
    energy = x*J_MPF_BML*x';   % engery of initial sequence 
end
Energy_trace = [energy];

while n<N
    xnew =x;
    energy_new = energy;
    pos = randi(stream,positions);
    index = randi(stream,[0,num_mutants_combine_array(1,pos)]); % index=0 i.e. consensus aa
    Start_site = num_mutants_combine_array_acc_all(pos)+1;
    End_site = num_mutants_combine_array_acc_all(pos+1);
    one = find(xnew(Start_site:End_site)>0);     % 目前seq x 在pos上mutant所在的位置
    if(~isempty(one) && ~logical(USE_H)) % 如果x在pos上不是consensus aa
        energy_new = energy_new-2*xnew*J_MPF_BML(:,Start_site+one-1)+J_MPF_BML(Start_site+one-1,Start_site+one-1);
    end
    if(index==0) 
        if(USE_H)
            energy_new = energy_new-xnew(Start_site:End_site)*H(Start_site:End_site)';
            xnew(Start_site:End_site)=0;                        
        else            
            xnew(Start_site:End_site)=0; 
        end
              
    else
        if(USE_H)
            energy_new = energy_new-xnew(Start_site:End_site)*H(Start_site:End_site)';
        end
        xnew(Start_site:End_site)=0;
        xnew(Start_site+index-1)=1;
         if(USE_H)
            energy_new = energy_new+xnew(Start_site+index-1)*H(Start_site+index-1)';
         else
            energy_new = energy_new+2*xnew*J_MPF_BML(:,Start_site+index-1)-J_MPF_BML(Start_site+index-1,Start_site+index-1); 
         end
    end

    trans_rate = exp(energy-energy_new);
    %trans_rate = 1/(1+exp((-energy+energy_new)));
    u = rand(stream);
    if energy_new<=energy
        x=xnew;
        energy=energy_new;
    elseif u < trans_rate
        x=xnew;
        energy=energy_new;
    end

    % sample
    %if(n+1>Burn_in)
        %if(mod((n+1),Thin_in)==0)
            %samples=[samples;x];
        %end
    %end
    Energy_trace = [Energy_trace; energy];
    n=n+1;
%     if rand < 1/(1+exp(xnew*triu(J_MPF_BML)*xnew'-x*triu(J_MPF_BML)*x'))

end
toc;
%% plots

% Energy plot
plot(0:N, Energy_trace)
xtitle('iteration')
ytitle('Sequence Energy')
title('MCMC for E2 1a sequence')

%% MCMC Metroplis_ revised by Haiyi
rng(0)
stream = RandStream("mcg16807","Seed",0);
tic;
samples=[];
N=1000; %1e7
Burn_in=0; %1e4
Thin_in=1; %1e3
USE_H=0;

num_mutants_combine_array_acc = cumsum(num_mutants_combine_array);
num_mutants_combine_array_acc_all = [0 num_mutants_combine_array_acc];
positions = size(num_mutants_combine_array,2);    
msa_bin_unique = unique(msa_bin,'row');
randvalue = randi(stream,[1 size(msa_bin_unique,1)]);
x = msa_bin_unique(randvalue,:);      % randomly pick a sequence from unique MSA sequence
n=0;

if(USE_H)
    energy = Calculate_Energy(x, J_MPF_BML,H);
else 
    energy = Calculate_Energy(x, J_MPF_BML);
end

while n<N
    xnew =x;
    pos = randi(stream,positions);
    index = randi(stream,[0,num_mutants_combine_array(1,pos)]);  % new mutant index=0 i.e. consensus aa
    Start_site = num_mutants_combine_array_acc_all(pos)+1;
    End_site = num_mutants_combine_array_acc_all(pos+1);
    xnew(Start_site:End_site)=0;           
    if index~=0
        xnew(Start_site+index-1)=1;
    end

    if(USE_H)
    energy_new = Calculate_Energy(xnew, J_MPF_BML,H);
    else 
    energy_new = Calculate_Energy(xnew, J_MPF_BML);
    end
   
    trans_rate = exp(energy-energy_new);
    %trans_rate = 1/(1+exp((-energy+energy_new)));
    u= rand(stream);
    if energy_new<=energy
        x=xnew;
        energy=energy_new;
    elseif u < trans_rate
        x=xnew;
        energy=energy_new;
    end

    % sample
    if(n+1>Burn_in)
        if(mod((n+1),Thin_in)==0)
            samples=[samples;x];
        end
    end

    n=n+1;
%     if rand < 1/(1+exp(xnew*triu(J_MPF_BML)*xnew'-x*triu(J_MPF_BML)*x'))

end
toc;
%% Single Mutation


num_mutants_combine_array_acc = cumsum(num_mutants_combine_array);


num_mutants_combine_array_acc_all = [0 num_mutants_combine_array_acc];
N = size(samples,1);
sites = size(msa_bin,2);
Single_Mutation_Model = sum(samples,1)/N;
Single_Mutation_Oberved = sum(msa_bin.*weight_seq,1)/sum(weight_seq);




figure;
plot(Single_Mutation_Oberved,Single_Mutation_Model,'*');
hold on;
arrayline_min = min(Single_Mutation_Oberved);
arrayline_max = max(Single_Mutation_Oberved);
arrayline = arrayline_min:0.01:arrayline_max;
plot(arrayline,arrayline,'k')



xlabel("P(i)-Observe");
ylabel("P(i)-Sample");




%% Double mutation


Double_Mutation_Oberved = zeros(sites,sites);
Double_Mutation_Model = Double_Mutation_Oberved;

for i =1:sites
    for j = i+1:sites
            num=0;
            for k =1:size(msa_bin,1)
                num=num+((msa_bin(k,i)==1)&&(msa_bin(k,j)==1))*weight_seq(k,1);           
            end
            Double_Mutation_Oberved(i,j)=num;
    end
end

Double_Mutation_Oberved = Double_Mutation_Oberved/sum(weight_seq);

Double_Mutation_Oberved_flat = Double_Mutation_Oberved(triu(true(size(Double_Mutation_Oberved)),1))';


for i =1:sites
    for j = i+1:sites
            num=0;
            for k =1:N
                num=num+((samples(k,i)==1)&&(samples(k,j)==1));           
            end
            Double_Mutation_Model(i,j)=num;
    end
end
Double_Mutation_Model = Double_Mutation_Model/N;

Double_Mutation_Model_flat = Double_Mutation_Model(triu(true(size(Double_Mutation_Model)),1))';



Index = true(sites,sites);
Index(1:num_mutants_combine_array_acc(1),1:num_mutants_combine_array_acc(1))=0;
for i =2:size(num_mutants_combine_array,2)
    Index(num_mutants_combine_array_acc(i-1)+1:num_mutants_combine_array_acc(i),num_mutants_combine_array_acc(i-1)+1:num_mutants_combine_array_acc(i))=0;
end

index=[];

for i =1:sites
   index = [index Index(i,i+1:sites)]; 
end

index = ~logical(index);
Double_Mutation_Oberved_flat(index)=[];
Double_Mutation_Model_flat(index)=[];


figure;
plot(Double_Mutation_Oberved_flat,Double_Mutation_Model_flat,'*')
hold on;
arrayline_min = min(Double_Mutation_Oberved_flat);
arrayline_max = max(Double_Mutation_Oberved_flat);
arrayline = arrayline_min:0.01:arrayline_max;
plot(arrayline,arrayline,'k')
xlabel("P(i,j)-Observe");
ylabel("P(i,j)-Sample");
%% Connect correlation
% Double_Mutation_Oberved = zeros(sites,sites);
% Double_Mutation_Model = Double_Mutation_Oberved;
% 
% for i =1:sites
%     for j = i+1:sites
%             num=0;
%             for k =1:size(msa_bin,1)
%                 num=num+((msa_bin(k,i)==1)&&(msa_bin(k,j)==1))*weight_seq(k,1);           
%             end
%             Double_Mutation_Oberved(i,j)=num;
%     end
% end
% Double_Mutation_Oberved = Double_Mutation_Oberved/sum(weight_seq);
% 
% 
% for i =1:sites
%     for j = i+1:sites
%             num=0;
%             for k =1:N
%                 num=num+((samples(k,i)==1)&&(samples(k,j)==1));           
%             end
%             Double_Mutation_Model(i,j)=num;
%     end
% end
% Double_Mutation_Model = Double_Mutation_Model/N;

% find location of two points in flattened array
zeros_remove = (ones(sites,sites));
 cumul_num_mutants_combine_array= cumsum(num_mutants_combine_array);
phi_extend = [1 cumul_num_mutants_combine_array+1];
for aba=1:length(phi_extend)-1
    array_sites = phi_extend(aba):phi_extend(aba+1)-1;
    for indi=1:length(array_sites)
        for indj=1:length(array_sites)
            if (indi~=indj)
                zeros_remove(array_sites(indi),array_sites(indj))=0;
            end
        end
    end
end
zeros_remove = sparse(zeros_remove);

ones_double = sparse(triu(ones(sites,sites),1));
double_pick_final = zeros_remove.*ones_double;
double_pick_final_flat = sparse(double_pick_final(:));
ind_double_pick_final = find(double_pick_final_flat==1); % location of two points

% Double_Mutation_Oberved  = (Double_Mutation_Oberved  +Double_Mutation_Oberved') - diag(diag(Double_Mutation_Oberved));
% Double_Mutation_Model  = (Double_Mutation_Model +Double_Mutation_Model') - diag(diag(Double_Mutation_Model));
connected_MCMC = Double_Mutation_Model - Single_Mutation_Model' *Single_Mutation_Model;
connected_MSA = Double_Mutation_Oberved  - Single_Mutation_Oberved'*Single_Mutation_Oberved;
connected_MCMC_flat = connected_MCMC(:);
connected_MSA_flat = connected_MSA(:);

connected_MCMC_flat = connected_MCMC_flat(ind_double_pick_final);
connected_MSA_flat = connected_MSA_flat(ind_double_pick_final);

figure;

arrayline_min = min(connected_MCMC_flat);
arrayline_max = max(connected_MCMC_flat);
arrayline = arrayline_min:0.01:arrayline_max;

plot(connected_MSA_flat,connected_MCMC_flat,'*');hold on;grid off;
plot(arrayline,arrayline,'k')
xlabel('Connected correlation (MSA)')
ylabel('Connected correlation')
%% Num of Mutations

Mutation_MSA = zeros(size(msa_aa,1),1);
Mutation_MCMC = zeros(size(msa_aa,1),1);

MSA=sum(msa_bin,2);
MC=sum(samples,2);

for i =1:size(msa_aa,1)
    for j =1:size(MSA,1)
        if(i==MSA(j))
            Mutation_MSA(i)=Mutation_MSA(i)+1;
        end
    end
end

Mutation_MSA = Mutation_MSA/size(msa_bin,1);

for i =1:size(msa_aa,1)
    for j =1:size(MC,1)
        if(i==MC(j))
            Mutation_MCMC(i)=Mutation_MCMC(i)+1;
        end
    end
end

Mutation_MCMC = Mutation_MCMC/size(samples,1);

left = min(min(MSA),min(MC));
right = max(max(MSA),max(MC));

Mutation_MCMC = Mutation_MCMC(left:right);
Mutation_MSA = Mutation_MSA(left:right);

figure;
semilogy(Mutation_MSA);
hold on;
semilogy(Mutation_MCMC);

legend('MSA','MCMC')

toc;
%% Calculate deltaE
num_mutants_combine_array_acc = cumsum(num_mutants_combine_array);
tic;

num_mutants_combine_array_acc_all = [0 num_mutants_combine_array_acc];
positions = size(num_mutants_combine_array,2);    
deltaE = zeros(1,size(samples,2));


P_samples = zeros(size(samples,1),1);

for i =1:size(samples,1)
    if(USE_H)
        P_samples(i) = samples(i,:)*H';
    else
         P_samples(i) = samples(i,:)*J_MPF_BML*samples(i,:)';
    end

end

P_samples = exp(-P_samples) ./sum(exp(-P_samples ));


for i =1:positions
    
    % wild type at that position
    ind_samples = find(all(samples(:,[num_mutants_combine_array_acc_all(i)+1:num_mutants_combine_array_acc_all(i+1)])==0,2));
    tmp_samples = samples(ind_samples,:);
    
    % sum over each sample
    for k =1:size(tmp_samples,1)
%         if(USE_H)           
%             E=tmp_samples(k,:)*H';
%         else
%             E=tmp_samples(k,:)*J_MPF_BML*tmp_samples(k,:)';
%         end
        for j =1:num_mutants_combine_array(i)
            index=num_mutants_combine_array_acc_all(i)+j;
            tmp_unique = tmp_samples(k,:);
            tmp_unique(index)=1;
            if(USE_H)
%                 deltaE(num_mutants_combine_array_acc_all(i)+j) = deltaE(num_mutants_combine_array_acc_all(i)+j)+(tmp_unique*H'-E)*weight/size(samples,1);
%                 deltaE(index) = deltaE(index)+H(index)/size(samples,1);
                deltaE(index) = deltaE(index)+H(index).*P_samples(ind_samples(k));
            else
%                 deltaE(num_mutants_combine_array_acc_all(i)+j) = deltaE(num_mutants_combine_array_acc_all(i)+j)+(tmp_unique*J_MPF_BML*tmp_unique'-E)/size(samples,1);
%                 deltaE(index) = deltaE(index)+(2*tmp_unique*J_MPF_BML(:,index)-J_MPF_BML(index,index))/size(samples,1);
                deltaE(index) = deltaE(index)+(2*tmp_unique*J_MPF_BML(:,index)-J_MPF_BML(index,index)).*P_samples(ind_samples(k));
            end


        end
    end
end

% %% Calculate average E

Avg_E = zeros(1,positions);
for i =1:positions
    delta = deltaE([num_mutants_combine_array_acc_all(i)+1:num_mutants_combine_array_acc_all(i+1)]);
    Avg_E(i) = sum(delta.*exp(-delta))/sum(exp(-delta));
    
end
toc;

%% Peak Analysis
Peaks = [];
tic;
num_mutants_combine_array_acc = cumsum(num_mutants_combine_array);
num_mutants_combine_array_acc_all = [0 num_mutants_combine_array_acc];
for i =1:size(msa_aa,1)
    
    x = msa_bin(i,:);

    x_E=x*J_MPF_BML*x';

    while 1
        lowest_E=x_E;
        for j =1:size(msa_aa,2)
            Start_site = num_mutants_combine_array_acc_all(j)+1;
            End_site   = num_mutants_combine_array_acc_all(j+1);
            one = find(x(Start_site:End_site)>0);
            if(one)
                zero_E = x_E-2*x*J_MPF_BML(:,Start_site+one-1)+J_MPF_BML(Start_site+one-1,Start_site+one-1);

            else
                zero_E = x_E;
            end
            low_E=zero_E;
            pos_start=Start_site;
            pos_end=End_site;
            pos_change=-1;
            for k = Start_site:End_site
                tmp_x=x;
                tmp_x(Start_site:End_site)=0;
                tmp_x(k)=1;
                E = zero_E+2*tmp_x*J_MPF_BML(:,k)-J_MPF_BML(k,k);
                if E<low_E
                    low_E = E;

                    pos_change=k;
                end
            end

            if low_E<lowest_E
                lowest_E = low_E;
                pos_start_all=pos_start;
                pos_end_all=pos_end;
                pos_change_all=pos_change;
            end
        end

        if abs(x_E-lowest_E)>1e-7
            x_E=lowest_E;
            if pos_change_all>0
               x(pos_start_all:pos_end_all)=0;
               x(pos_change_all)=1;
               
            else
               x(pos_start_all:pos_end_all)=0;
            end
        else
            Peaks=[Peaks;x];
            break;
        end
        
    end
    

end
toc;


%% Contact Prediction (Frobenius norm)
tic;

num_mutants_combine_array_acc = cumsum(num_mutants_combine_array);


num_mutants_combine_array_acc_all = [0 num_mutants_combine_array_acc];
num = size(msa_aa,2);
F = zeros(num,num);
% J_MPF_BML = 2*J_MPF_BML-diag(diag(J_MPF_BML));
for i = 1:num
    Start_site_i = num_mutants_combine_array_acc_all(i)+1;
    End_site_i   = num_mutants_combine_array_acc_all(i+1);
   
    for j =1:num
        Start_site_j = num_mutants_combine_array_acc_all(j)+1;
        End_site_j   = num_mutants_combine_array_acc_all(j+1); 
        F(i,j)=F(i,j)+sum(sum((J_MPF_BML(Start_site_i:End_site_i,Start_site_j:End_site_j).^2)));
    
    end
end
% F = (F+F')/2;
F = F.^0.5;
F_new=F;
avg=0;
for i = 1:num  
    for j =1:num
        if i ~=j
            avg = avg+F(i,j);
        end
    end
end

% avg = mean(mean(F));
avg = avg/(num^2-num);



for i =1:num
    for j =1:num
        
          avg_i = mean(F(:,j))-F(j,j)/num;
          avg_j = mean(F(i,:))-F(i,i)/num;
%         avg_i = mean(F(:,j));
%         avg_j = mean(F(i,:));
        F_new(i,j) = F(i,j) - avg_i*avg_j/avg;
    end
end
F = F_new;
F =triu(F,1);


Nmax = num*(num-1)/2; % get Nmax biggest entries
[ ~, Ind ] = sort(F(:),1,'descend');

[ ind_row, ind_col ] = ind2sub(size(F),Ind(1:Nmax)); % fetch indices
load('6mei.mat')
indx=indx-383;
indy=indy-383;

% load('I-TASSER\model5.mat')
ind_non_conserve = setdiff(1:363,conserved);
ind_col = ind_non_conserve(ind_col);
ind_row = ind_non_conserve(ind_row);

ia1 = ind_row<=(647-383) & ind_row>= (410-383);
ia2 = ind_col<=(647-383) & ind_col>= (410-383);
ia = ia1 & ia2;
ind_col =ind_col(ia);
ind_row =ind_row(ia);
% Ind = or(ind_col<ind_row,ind_col>(max(indy)-sum(conserved<=max(indy))) ); 
% % Ind = abs(ind_col-ind_row)<=5|ind_col<ind_row|ind_col>(max(indy)-sum(conserved<=max(indy)))  ;
% 
% ind_col(Ind)=[];
% ind_row(Ind)=[];

% [~,ia,~]=intersect(indx,conserved);
% indx(ia)=[];
% indy(ia)=[];
% for i = 1:length(indx)
%     indx(i) = indx(i)-sum(conserved<=indx(i));
%     indy(i) = indy(i)-sum(conserved<=indy(i));
% end


TPR = zeros(size(ind_row));



count =0;
for i =1:length(ind_row)
   ia1 = find(indx==ind_row(i));
   ia2 = find(indy==ind_col(i));
   if(intersect(ia1,ia2))
        count = count+1;
   end

    TPR(i) = count/i;

end




FIG=figure;
plot(TPR,'Color','k','LineWidth',1)
hold on


set(gca, 'XScale', 'log')
xlim([1 1e3])

xlabel('Top x pairs')
ylabel({'Precision'})
FIG.Units = 'centimeters';
FIG.Name = 'Ours';


FIG.Units = 'centimeters';
% set(gcf,'Position',[10 10 7.84 6]);
% set(gca,'Position',[.15 .2 .8 .74]);  %调整 XLABLE和YLABLE不会被切掉
set(gcf,'Position',[10 10 7.84 6]);
set(gca,'Position',[.15 .2 .8 .74]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600');
%% DCA
tic;
load('DCA_Ahmed.mat');
num=363;
F = DI_new;
F =triu(F,1);
Nmax = num*(num-1)/2; % get Nmax biggest entries
[ ~, Ind ] = sort(F(:),1,'descend');
[ ind_row, ind_col ] = ind2sub(size(F),Ind(1:Nmax)); % fetch indices
load('6mei.mat')
indx=indx-383;
indy=indy-383;
ia1 = ind_row<=(647-383) & ind_row>= (410-383);
ia2 = ind_col<=(647-383) & ind_col>= (410-383);
ia = ia1 & ia2;
ind_col =ind_col(ia);
ind_row =ind_row(ia);


% load('I-TASSER\model5.mat')
TPR = zeros(size(ind_row));
count =0;
for i =1:length(ind_row)
   ia1 = find(indx==ind_row(i));
   ia2 = find(indy==ind_col(i));
   if(intersect(ia1,ia2))
        count = count+1;
   end

    TPR(i) = count/i;

end

FIG=figure
color = [0.6000    0.6000    0.6000];%gray
plot(TPR,'Color','k','LineWidth',1);

% 
% box off
set(gca, 'XScale', 'log')
xlim([1 1e3])

xlabel('Top x pairs')
ylabel({'Precision'})
FIG.Units = 'centimeters';
FIG.Name = 'DCA';

FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 10.45 7.84]);
set(gca,'Position',[.15 .2 .8 .74]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600');

%%
load('E2.mat')
ind_col =indy;
ind_row =indx;

load('6mei.mat')
indx=indx-383;
indy=indy-383;
ia1 = ind_row<=(647-383) & ind_row>= (410-383);
ia2 = ind_col<=(647-383) & ind_col>= (410-383);
ia = ia1 & ia2;
ind_col =ind_col(ia);
ind_row =ind_row(ia);

% load('I-TASSER\model5.mat')
TPR = zeros(size(ind_row));
count =0;
for i =1:length(ind_row)
   ia1 = find(indx==ind_row(i));
   ia2 = find(indy==ind_col(i));
   if(intersect(ia1,ia2))
        count = count+1;
   end

    TPR(i) = count/i;

end

FIG=figure
color = [0.6000    0.6000    0.6000];%gray
plot(TPR,'Color','k','LineWidth',1);

% 
% box off
set(gca, 'XScale', 'log')
xlim([1 1e3])

xlabel('Top x pairs')
ylabel({'Precision'})
FIG.Units = 'centimeters';
FIG.Name = 'DCA';

FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 10.45 7.84]);
set(gca,'Position',[.15 .2 .8 .74]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600');