load('results/workspace/test_1226.mat')

Single_Mutation_Observed = sum(msa_bin.*weight_seq,1)/sum(weight_seq);  % f_i(a)
H = -log((Single_Mutation_Observed)./(1-Single_Mutation_Observed));
conserve=0;
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

Energy1_norm = Energy1 - Energy1(1);
Energy2_norm = Energy2 - Energy2(1);
Energy3_norm = Energy3 - Energy3(1);

fitness1_norm = fitness1 ./ fitness1(1);
fitness2_norm = fitness2 ./ fitness2(1);
fitness3_norm = fitness3 ./ fitness3(1);

%%  Prentoe2019 group4
[  ~ , test1_seq ] = fastaread('data/Prentoe2019.fasta');
test1_seq = cell2mat(test1_seq');
test1_seq_bin = Binary_Seq(test1_seq,amino_single_combine_array,ind_conserve);

if conserve == 1
    Energy4 = Calculate_Energy(test1_seq_bin, J_MPF_BML,H);
else 
    Energy4 = Calculate_Energy(test1_seq_bin, J_MPF_BML);
end

Energy4_norm = Energy4 - Energy4(1);

fitness4 = [100, 144.9275, 129.4686, 8.696, 43.478, 133.3333, 10.628, 38.6473, 84.058]' ;
fitness4_norm = fitness4 ./ fitness4(1);
%% 
markersize = 6;
line_width = 0.5;

all_fitness = [fitness1_norm' fitness2_norm' fitness3_norm' fitness4_norm' ];
all_Energy = [Energy1_norm', Energy2_norm', Energy3_norm' Energy4_norm' ];

P = polyfit (all_Energy, all_fitness, 1);
x = -5 :0.5:20;
y = P(1)*x +P(2);

figure
hold on 
plot(all_Energy,all_fitness,'o','MarkerSize',markersize)
plot(x,y,'k--','LineWidth',1);
text(-2.5 ,-0.25,sprintf('$$r= %.4f $$', P(1)), 'Interpreter','latex', 'FontSize', 16)
P(1)
xlabel('Energy (relatived)')
ylabel('Experimental fitness (relatived)')
%legend('Esteban-Riesco2013-1,Esteban-Riesco2013-2,Esteban-Riesco2013-3,Prentoe2019')
title('calculation of b')
hold off




