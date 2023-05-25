function []= Biological_validation_3a(msa_bin,weight_seq,amino_single_combine_array,ind_conserve,J_MPF_BML,conserve)
%%
markersize = 6;
line_width = 0.5;

% conserved-only model
Single_Mutation_Observed = sum(msa_bin.*weight_seq,1)/sum(weight_seq);  % f_i(a)
H = -log((Single_Mutation_Observed)./(1-Single_Mutation_Observed));

%% Riesco2013 group1-3

[  ~ , test1_seq ] = fastaread('data/Riesco2013.fasta');
test1_seq = cell2mat(test1_seq');
test1_seq_bin = Binary_Seq(test1_seq,amino_single_combine_array,ind_conserve);

if conserve == 1
    Energy1 = Calculate_Energy(test1_seq_bin(1:2,:), J_MPF_BML,H);
else 
    Energy1 = Calculate_Energy(test1_seq_bin(1:2,:), J_MPF_BML);
end

if conserve == 1
    Energy2 = Calculate_Energy(test1_seq_bin(3:5,:), J_MPF_BML,H);
else 
    Energy2 = Calculate_Energy(test1_seq_bin(3:5,:), J_MPF_BML);
end

if conserve == 1
    Energy3 = Calculate_Energy(test1_seq_bin(6:7,:), J_MPF_BML,H);
else 
    Energy3 = Calculate_Energy(test1_seq_bin(6:7,:), J_MPF_BML);
end

fitness1 = [52117.21, 8404.818]';
fitness2 = [4015.814, 10907.77, 95.36334]';
fitness3 = [199.0532, 95.36334]';

Energy1_norm = normalize(Energy1);
Energy2_norm = normalize(Energy2);
Energy3_norm = normalize(Energy3);

fitness1_norm = normalize(fitness1);
fitness2_norm = normalize(fitness2);
fitness3_norm = normalize(fitness3);

[r1,p1] = corr(Energy1_norm , fitness1_norm , 'type', 'Spearman');
[r2,p2] = corr(Energy2_norm , fitness2_norm , 'type', 'Spearman');
[r3,p3] = corr(Energy3_norm , fitness3_norm , 'type', 'Spearman');

if conserve == 1
    title('Conservational-only')
else 
    title('Our Model')
end

%% Prentoe2019 group4
[  ~ , test1_seq ] = fastaread('data/Prentoe2019.fasta');
test1_seq = cell2mat(test1_seq');
test1_seq_bin = Binary_Seq(test1_seq,amino_single_combine_array,ind_conserve);

if conserve == 1
    Energy4 = Calculate_Energy(test1_seq_bin, J_MPF_BML,H);
else 
    Energy4 = Calculate_Energy(test1_seq_bin, J_MPF_BML);
end

Energy4_approx=Energy4(1:9);
Energy4_norm = normalize(Energy4_approx);

fitness4 = [100, 144.9275, 129.4686, 8.696, 43.478, 133.3333, 10.628, 38.6473, 84.058]' ;
fitness4_approx =fitness4(1:9);
fitness4_norm = normalize(fitness4_approx);

[r4,p4] = corr(Energy4_norm , fitness4_norm , 'type', 'Spearman');

if conserve == 1
    title('Conservational-only')
else 
    title('Our Model')
end

%%
all_fitness = [fitness1_norm' fitness2_norm' fitness3_norm' fitness4_norm' ];
all_Energy = [Energy1_norm', Energy2_norm', Energy3_norm' Energy4_norm' ];

total_length = length( all_fitness);
w1 = length(fitness1') / total_length;
w2 = length(fitness2') / total_length;
w3 = length(fitness3') / total_length;
w4 =  length(fitness4_approx') / total_length;

ave_r = [w1 w2 w3 w4] * [r1 r2 r3 r4]';

P = polyfit (all_Energy, all_fitness, 1);
x = -3 :0.5:3;
y = P(1)*x +P(2);

figure
hold on 
plot(Energy1_norm,fitness1_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','cyan', 'LineWidth',line_width)
plot(Energy2_norm,fitness2_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','green', 'LineWidth',line_width)
plot(Energy3_norm,fitness3_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','blue', 'LineWidth',line_width)
plot(Energy4_norm , fitness4_norm, 'o' , 'MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','red', 'LineWidth',line_width)
plot(x,y,'k--','LineWidth',1);
xlabel('Energy (normalized)')
ylabel('Experimental fitness (normalized)')
text(-1.75 ,-1,sprintf('$$r= %.4f $$', ave_r), 'Interpreter','latex', 'FontSize', 16)
legend('Esteban-Riesco2013-1', 'Esteban-Riesco2013-2', 'Esteban-Riesco2013-3', 'Prentoe2019')

if conserve == 1
    title('Conservational-only')
else 
    title('Our Model')
end
hold off

end