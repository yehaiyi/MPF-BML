clear;
load('results/workspace/MPF_BML_output.mat')
%run startup.m

%% Serre2013
[  ~ , test1_seq ] = fastaread('data/Serre2013.fasta');
test1_seq = cell2mat(test1_seq');
test1_seq_bin = Binary_Seq(test1_seq,amino_single_combine_array,ind_conserve);

Energy1 = Calculate_Energy(test1_seq_bin, J_MPF_BML);
fitness1 = [132629.9205, 183717.6762]';

Energy1_norm = normalize(Energy1);
fitness1_norm = normalize(fitness1);

[r1,p1] = corr(Energy1_norm , fitness1_norm , 'type', 'Spearman')

markersize = 6;
line_width = 0.5;
figure;
plot(Energy1_norm,fitness1_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','red', 'LineWidth',line_width)

xlabel('Energy (normalized)')
ylabel('Experimental fitness (normalized)')

%% Urbanowicz2016
[  ~ , test1_seq ] = fastaread('data/Urbanowicz2016.fasta');
test1_seq = cell2mat(test1_seq');
test1_seq_bin = Binary_Seq(test1_seq,amino_single_combine_array,ind_conserve);

Energy1 = Calculate_Energy(test1_seq_bin, J_MPF_BML);
fitness1 = [628.7927, 2022.88]';

Energy1_norm = normalize(Energy1);
fitness1_norm = normalize(fitness1);

[r1,p1] = corr(Energy1_norm , fitness1_norm , 'type', 'Spearman')

markersize = 6;
line_width = 0.5;
figure;
plot(Energy1_norm,fitness1_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','red', 'LineWidth',line_width)

xlabel('Energy (normalized)')
ylabel('Experimental fitness (normalized)')

%% Riesco2013
markersize = 6;
line_width = 0.5;

[  ~ , test1_seq ] = fastaread('data/Riesco2013.fasta');
test1_seq = cell2mat(test1_seq');
test1_seq_bin = Binary_Seq(test1_seq,amino_single_combine_array,ind_conserve);

Energy1 = Calculate_Energy(test1_seq_bin(1:2,:), J_MPF_BML);
Energy2 = Calculate_Energy(test1_seq_bin(3:5,:), J_MPF_BML);
Energy3 = Calculate_Energy(test1_seq_bin(6:7,:), J_MPF_BML);

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

figure
hold on 
plot(Energy1_norm,fitness1_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','red', 'LineWidth',line_width)
P1 = polyfit(Energy1_norm, fitness1_norm, 1 );
x1 = -1:0.5:1.5;
y1 = P1(1)*x1 +P1(2);
plot(x1,y1,'k--','LineWidth',1);
xlabel('Energy (normalized)')
ylabel('Experimental fitness (normalized)')
text(-0.75,-1,sprintf('$$r= %.2f $$', r1), 'Interpreter','latex', 'FontSize', 16)
legend('Esteban-Riesco2013-1')
hold off

figure
hold on 
plot(Energy2_norm,fitness2_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','green', 'LineWidth',line_width)
P2 = polyfit(Energy2_norm, fitness2_norm, 1 );
x2 = -1:0.5:1.5;
y2 = P2(1)*x2 +P2(2);
plot(x2,y2,'k--','LineWidth',1);
xlabel('Energy (normalized)')
ylabel('Experimental fitness (normalized)')
text(-0.75,-1,sprintf('$$r= %.2f $$', r2), 'Interpreter','latex', 'FontSize', 16)
legend('Esteban-Riesco2013-2')
hold off

figure
hold on 
plot(Energy3_norm,fitness3_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','blue', 'LineWidth',line_width)
P3 = polyfit(Energy3_norm, fitness3_norm, 1 );
x3 = -1:0.5:1.5;
y3 = P3(1)*x3 +P3(2);
plot(x3,y3,'k--','LineWidth',1);
xlabel('Energy (normalized)')
ylabel('Experimental fitness (normalized)')
text(-0.75,-1,sprintf('$$r= %.2f $$', r3), 'Interpreter','latex', 'FontSize', 16)
legend('Esteban-Riesco2013-3')
hold off

%% 
all_fitness = [fitness1_norm' fitness2_norm' fitness3_norm'];
all_Energy = [Energy1_norm', Energy2_norm', Energy3_norm'];
total_length = length( all_fitness);
w1 = length(fitness1') / total_length;
w2 = length(fitness2') / total_length;
w3 = length(fitness3') / total_length;

ave_r = [w1 w2 w3] * [r1 r2 r3]';

P = polyfit (all_Energy, all_fitness, 1);
x = -1:0.5:1.5;
y = P(1)*x +P(2);

figure
hold on 
plot(Energy1_norm,fitness1_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','red', 'LineWidth',line_width)
plot(Energy2_norm,fitness2_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','green', 'LineWidth',line_width)
plot(Energy3_norm,fitness3_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','blue', 'LineWidth',line_width)
plot(x,y,'k--','LineWidth',1);
xlabel('Energy (normalized)')
ylabel('Experimental fitness (normalized)')
text(-0.75,-1,sprintf('$$r= %.4f $$', ave_r), 'Interpreter','latex', 'FontSize', 16)
legend('Esteban-Riesco2013-1', 'Esteban-Riesco2013-2', 'Esteban-Riesco2013-3')








