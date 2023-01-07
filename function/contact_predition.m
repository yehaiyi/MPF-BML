%% load model
load('results/workspace/1215.mat')
load data\weight_seq.mat
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

%% top x pairs 对应的precision
load data\contact_mat.mat


top_1000_contact = Top_x_pairs_contact(F_APC,N,1000,contact_mat);
top_1000_precision = zeros(1,1000);

for i=1:1000
    top_1000_precision(i) = sum(top_1000_contact(1:i)) / i;
end

hold on 
plot(log10(1:10),top_1000_precision(1:10),'Marker','.','MarkerSize',10,'Color','k')
plot(log10(1:1000),top_1000_precision)
xticks(log10(1:10))
xticks(log10([1 10 100 1000]))
xticklabels({'10^{0}','10^{1}','10^{2}','10^{3}'})
xlabel('Top x pairs')
ylabel('Precision')
legend("Precision of top 10 pairs")
title("Precision of contact predictions on 3a E2 residues")





