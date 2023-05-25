%% load model
load('results/workspace/test_1226.mat')  %3a model
%% 

load('data/Model_1a.mat')
%load data\weight_seq.mat
%% 
load('data/Model_1b.mat')
ind_non_conserve = find (mutantfreq ~= 0)
%% 计算F_APC
N = length(ind_non_conserve);
Frob_J = zeros(N,N);
F_APC = zeros(N,N)-10;


for i= 1: N
    for j= 1: N
        [i_s, i_e] = res2col(amino_single_combine_array,i);
        [j_s, j_e] = res2col(amino_single_combine_array,j);
        J_ij = J_MPF_BML(i_s:i_e, j_s:j_e);
        Frob_J(i,j) = norm(J_ij,'fro');
    end
end


for i= 1: N
    for j= 1: N
        if i>= j
            continue
        end
        F_i = (sum(Frob_J(i,:))-Frob_J(i,i) ) / (N-1) ;
        F_j = (sum(Frob_J(:,j))-Frob_J(j,j) ) / (N-1) ;
        F = ( sum(Frob_J,'all') - sum(diag(Frob_J)) )/ (N^2-N) ; 
        F_APC(i,j) = Frob_J(i,j) - F_i*F_j / F;
    end
end

%% top x pairs 对应的precision 针对contact mat是完整的
load data\H77E2_contact_mat.mat

top_1000_contact_H77 = Top_x_pairs_contact(F_APC,N,1000,H77E2_contact_mat,ind_non_conserve);
top_1000_precision_H77 = zeros(1,1000);

for i=1:1000
    top_1000_precision_H77(i) = sum(top_1000_contact_H77(1:i)) / i;
end

figure
hold on 
plot(log10(1:1000),top_1000_precision_H77,'Color',blue)
plot(log10(1:10),top_1000_precision_H77(1:10),'Marker','.','MarkerSize',10,'Color','k')

xticks(log10(1:10))
xticks(log10([1 10 100 1000]))
xticklabels({'10^{0}','10^{1}','10^{2}','10^{3}'})
legend('H77 Alphafold protein')
xlabel('Top x pairs')
title("Precision of contact predictions on 1a E2(with all the residues)")

%% 当contact mat不完整，e.g. 4mwf， 只考虑有contact信息的那些residues
% pairs是固定的，只是不同contact map的区别
load data\4mwf_contact_mat_8.mat
load data\4mwf_real_index.mat real_index
load data\H77E2_contact_mat.mat

[top_1000_contact_4mwf, top_value_real_rank] = Top_x_pairs_contact_incompleted(F_APC,N,1000,contact_mat,ind_non_conserve,real_index);
[top_1000_contact_H77p, ~] = Top_x_pairs_contact_incompleted(F_APC,N,1000,H77E2_contact_mat,ind_non_conserve,real_index);
%% 
top_1000_precision_4mwf = zeros(1,1000);
top_1000_precision_H77p = zeros(1,1000);


for i=1:1000
    top_1000_precision_4mwf(i) = sum(top_1000_contact_4mwf(1:i)) / i;
    top_1000_precision_H77p(i) = sum(top_1000_contact_H77p(1:i)) / i;
end

figure
hold on 
plot(log10(1:1000),top_1000_precision_4mwf,'Color',orange)
plot(log10(1:1000),top_1000_precision_H77p,'Color',blue)
plot(log10(1:10),top_1000_precision_4mwf(1:10),'Marker','.','MarkerSize',10)
plot(log10(1:10),top_1000_precision_H77p(1:10),'Marker','.','MarkerSize',10)


xticks(log10(1:10))
xticks(log10([1 10 100 1000]))
xticklabels({'10^{0}','10^{1}','10^{2}','10^{3}'})
xlabel('Top x pairs (within part of the residues)')
ylabel('Precision')
legend("Real protein(4mwf)",'H77 Alphafold protein')
title("Precision of contact predictions on 1a E2 (part of the residues)")

%% 6MEI 不完整contact mat
load data\6mei_contact_mat.mat
load data\6mei_real_index.mat
load data\AF_6mei_contact_mat.mat
load data\Con1_contact_mat.mat

x=500;

[top_1000_contact_6mei, top_value_real_rank] = Top_x_pairs_contact_incompleted(F_APC,N,x,contact_mat_6mei,ind_non_conserve,real_index_6mei);
%% 

[top_1000_contact_Con1p, ~] = Top_x_pairs_contact_incompleted(F_APC,N,x,Con1_contact_mat,ind_non_conserve,real_index_6mei);
%% 
top_1000_precision_6mei = zeros(1,x);
top_1000_precision_Con1p = zeros(1,x);


for i=1:x
    top_1000_precision_6mei(i) = sum(top_1000_contact_6mei(1:i)) / i;
    top_1000_precision_Con1p(i) = sum(top_1000_contact_Con1p(1:i)) / i;
end

figure
hold on 
plot(log10(1:x),top_1000_precision_6mei,'Color',orange)
plot(log10(1:x),top_1000_precision_Con1p,'Color',blue)
plot(log10(1:10),top_1000_precision_6mei(1:10),'Marker','.','MarkerSize',10)
plot(log10(1:10),top_1000_precision_Con1p(1:10),'Marker','.','MarkerSize',10)


xticks(log10(1:10))
xticks(log10([1 10 100 1000]))
xticklabels({'10^{0}','10^{1}','10^{2}','10^{3}'})
xlabel('Top x pairs (within part of the residues)')
ylabel('Precision')
legend("Real protein(6mei)",'Con1 Alphafold protein')
title("Precision of contact predictions on 1b E2 (part of the residues)")

%% 1b 完整residues的
load data\Con1_contact_mat.mat

top_1000_contact_Con1 = Top_x_pairs_contact(F_APC,N,1000,Con1_contact_mat,ind_non_conserve);
top_1000_precision_Con1 = zeros(1,1000);

for i=1:1000
    top_1000_precision_Con1(i) = sum(top_1000_contact_Con1(1:i)) / i;
end

figure
hold on 
plot(log10(1:1000),top_1000_precision_Con1,'Color',blue)
plot(log10(1:10),top_1000_precision_Con1(1:10),'Marker','.','MarkerSize',10,'Color','k')

xticks(log10(1:10))
xticks(log10([1 10 100 1000]))
xticklabels({'10^{0}','10^{1}','10^{2}','10^{3}'})
legend('Con1 Alphafold protein')
xlabel('Top x pairs')
ylabel('Precision')
title("Precision of contact predictions on 1b E2(with all the residues)")

%% random pairs
X = round(rand(1,1000)*339)+1;
Y  = round(rand(1,1000)*339)+1;

top_1000_contact = zeros(1,1000);

for i = 1:1000
        r_i = ind_non_conserve(X(i));
        c_i = ind_non_conserve(Y(i));
        top_1000_contact(i) = contact_mat(r_i,c_i);
end

random_1000_precision = zeros(1,1000);

for i=1:1000
    random_1000_precision(i) = sum(top_1000_contact(1:i)) / i;
end

hold on 
plot(log10(1:10),random_1000_precision(1:10),'Marker','.','MarkerSize',10,'Color','k')
plot(log10(1:1000),random_1000_precision)
xticks(log10(1:10))
xticks(log10([1 10 100 1000]))
xticklabels({'10^{0}','10^{1}','10^{2}','10^{3}'})
xlabel('Top x pairs')
ylabel('Precision')
legend("Precision of top 10 pairs")
title("Precision of contact predictions on 3a E2 residues")


%% 
load data\contact_mat.mat

top_1000_contact_S52 = Top_x_pairs_contact(F_APC,N,1000,contact_mat,ind_non_conserve);
top_1000_precision_S52 = zeros(1,1000);

for i=1:1000
    top_1000_precision_S52(i) = sum(top_1000_contact_S52(1:i)) / i;
end

figure
hold on 

plot(log10(1:1000),top_1000_precision_H77)
plot(log10(1:1000),top_1000_precision_Con1)
plot(log10(1:1000),top_1000_precision_S52)
%plot(log10(1:10),top_1000_precision_H77(1:10),'Marker','.','MarkerSize',10)
%plot(log10(1:10),top_1000_precision_Con1(1:10),'Marker','.','MarkerSize',10)
%plot(log10(1:10),top_1000_precision_S52(1:10),'Marker','.','MarkerSize',10)

xticks(log10(1:10))
xticks(log10([1 10 100 1000]))
xticklabels({'10^{0}','10^{1}','10^{2}','10^{3}'})
legend('1a (AF H77)','1b (AF Con1)','3a (AF S52)')
xlabel('Top x pairs')
ylabel('Precision')
title("Precision of contact predictions with 3 subtypes (AF:AlphaFold)")









