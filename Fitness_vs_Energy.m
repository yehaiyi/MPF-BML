load('results/workspace/test_1226.mat')
load data\weight_seq.mat
%run startup.m
%%
markersize = 6;
line_width = 0.5;

% conserved-only model
Single_Mutation_Observed = sum(msa_bin.*weight_seq,1)/sum(weight_seq);  % f_i(a)
H = -log((Single_Mutation_Observed)./(1-Single_Mutation_Observed));

conserve=0;

%% 导入369的msa_aa
load data\NumofPatient_3aE2.mat;
inputfile = 'data/3a_E2_ori.fasta';
[Header_fasta, Sequence_fasta] = fastaread(inputfile);
msa_aa_369 = cell2mat(Sequence_fasta');

% preprocess 
no_patient_idx = find(~patient);
if length(no_patient_idx) >0
    msa_aa_369(no_patient_idx) = []
end

load data\outliers.mat
msa_aa_369(outliers,:) =[];
patient(outliers,:) = [];


%% Serre2013
[  ~ , test1_seq ] = fastaread('data/Serre2013.fasta');
test1_seq = cell2mat(test1_seq');
test1_seq_bin = Binary_Seq(test1_seq,amino_single_combine_array,ind_conserve);


if conserve == 1
    Energy6 = Calculate_Energy(test1_seq_bin, J_MPF_BML,H);
else 
    Energy6 = Calculate_Energy(test1_seq_bin, J_MPF_BML);
end

fitness6 = [63506, 93719]';

Energy6_norm = normalize(Energy6);
fitness6_norm = normalize(fitness6);

[r6,p6] = corr(Energy6_norm , fitness6_norm , 'type', 'Spearman');

figure;
hold on
plot(Energy6_norm,fitness6_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','red', 'LineWidth',line_width)
P1 = polyfit(Energy6_norm, fitness6_norm, 1 );
x1 = -1:0.5:1.5;
y1 = P1(1)*x1 +P1(2);
plot(x1,y1,'k--','LineWidth',1);
xlabel('Energy (normalized)')
ylabel('Experimental fitness (normalized)')
text(-0.75,-1,sprintf('$$r= %.2f $$', r6), 'Interpreter','latex', 'FontSize', 16)
legend(' Serre2013')
hold off

if conserve == 1
    title('Conservational-only')
else 
    title('Our Model')
end

%% Urbanowicz2016

[  ~ , test1_seq ] = fastaread('data/Urbanowicz2016.fasta');
test1_seq = cell2mat(test1_seq');
test1_seq_bin = Binary_Seq(test1_seq,amino_single_combine_array,ind_conserve);

if conserve == 1
    Energy7 = Calculate_Energy(test1_seq_bin, J_MPF_BML,H);
else 
    Energy7 = Calculate_Energy(test1_seq_bin, J_MPF_BML);
end

fitness7 = [628.7927, 2022.88]';

Energy7_norm = normalize(Energy7);
fitness7_norm = normalize(fitness7);

[r7,p7] = corr(Energy7_norm , fitness7_norm , 'type', 'Spearman')

markersize = 6;
line_width = 0.5;
figure;
hold on
plot(Energy7_norm,fitness7_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','red', 'LineWidth',line_width)
P1 = polyfit(Energy7_norm, fitness7_norm, 1 );
x1 = -1:0.5:1.5;
y1 = P1(1)*x1 +P1(2);
plot(x1,y1,'k--','LineWidth',1);

xlabel('Energy (normalized)')
ylabel('Experimental fitness (normalized)')
legend(' Urbanowicz2016')

if conserve == 1
    title('Conservational-only')
else 
    title('Our Model')
end

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

figure
hold on 
plot(Energy1_norm,fitness1_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','cyan', 'LineWidth',line_width)
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
%Energy4_approx = [Energy4(1:2)'  Energy4(3) Energy4(5)  Energy4(6)  Energy4(8) ]' ;
%Energy4_approx = [Energy4(1:2)'  Energy4(3) Energy4(4) Energy4(5) Energy4(6) Energy4(7) Energy4(8) ]' ;
Energy4_norm = normalize(Energy4_approx);


fitness4 = [100, 144.9275, 129.4686, 8.696, 43.478, 133.3333, 10.628, 38.6473, 84.058]' ;
fitness4_approx =fitness4(1:9);
%fitness4_approx = [fitness4(1:2)'  fitness4(3) fitness4(5)  fitness4(6)  fitness4(8)]' ; 
%fitness4_approx = [fitness4(1:2)'  fitness4(3) fitness4(4)  fitness4(6)  fitness4(8) fitness4(9)]' ; 
fitness4_norm = normalize(fitness4_approx);

[r4,p4] = corr(Energy4_norm , fitness4_norm , 'type', 'Spearman');

figure
hold on
plot(Energy4_norm , fitness4_norm, 'o' , 'MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','red', 'LineWidth',line_width)
P4 = polyfit(Energy4_norm, fitness4_norm, 1 );
x4 = -3:0.5:3;
y4 = P4(1)*x4 +P4(2);
plot(x4,y4,'k--','LineWidth',1);
xlabel('Energy (normalized)')
ylabel('Experimental fitness (normalized)')
text(-1.75,-1,sprintf('$$r= %.2f $$', r4), 'Interpreter','latex', 'FontSize', 16)

%text(Energy4_norm(1),fitness4_norm(1)+0.1,'S52/JFH1')
%text(Energy4_norm(2),fitness4_norm(2)+0.1,'N417D')
%text(Energy4_norm(3),fitness4_norm(3)+0.1,'N430D')
%text(Energy4_norm(4),fitness4_norm(4)+0.1,'N448D')
%text(Energy4_norm(5),fitness4_norm(5)+0.1,'N556D')
%text(Energy4_norm(6),fitness4_norm(6)+0.1,'N623D')

if conserve == 1
    title('Conservational-only')
else 
    title('Our Model')
end




legend('Prentoe2019')
hold off


%%  Alzua2022 group5
[  ~ , test1_seq ] = fastaread('data/Alzua2022.fasta');
test1_seq = cell2mat(test1_seq');
test1_seq_bin = Binary_Seq(test1_seq,amino_single_combine_array,ind_conserve);

if conserve == 1
    Energy5 = Calculate_Energy(test1_seq_bin, J_MPF_BML,H);
else 
    Energy5 = Calculate_Energy(test1_seq_bin, J_MPF_BML);
end

fitness5 = [2935.905, 10378.37, 15615.23]';
%fitness5 = [2935.905, 10378.37, 15615.23]';
%fitness5 = [2935.905, 15615.23]';

Energy5_norm = normalize(Energy5);
fitness5_norm = normalize(fitness5);

[r5,p5] = corr(Energy5_norm , fitness5_norm , 'type', 'Spearman');

figure;
hold on
plot(Energy5_norm,fitness5_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','m', 'LineWidth',line_width)
P1 = polyfit(Energy5_norm, fitness5_norm, 1 );
x1 = -1.5:0.5:2;
y1 = P1(1)*x1 +P1(2);
plot(x1,y1,'k--','LineWidth',1);
xlabel('Energy (normalized)')
ylabel('Experimental fitness (normalized)')
text(-0.75,-1,sprintf('$$r= %.2f $$', r5), 'Interpreter','latex', 'FontSize', 16)
legend(' Alzua2022')
hold off

if conserve == 1
    title('Conservational-only')
else 
    title('Our Model')
end





%% 
all_fitness = [fitness1_norm' fitness2_norm' fitness3_norm' fitness4_norm' ];
%all_fitness = [fitness1_norm' fitness2_norm' fitness3_norm' fitness4_norm' fitness5_norm' fitness6_norm'];


all_Energy = [Energy1_norm', Energy2_norm', Energy3_norm' Energy4_norm' ];
%all_Energy = [Energy1_norm', Energy2_norm', Energy3_norm' Energy4_norm' Energy5_norm' Energy6_norm'];

total_length = length( all_fitness);
w1 = length(fitness1') / total_length;
w2 = length(fitness2') / total_length;
w3 = length(fitness3') / total_length;
w4 =  length(fitness4_approx') / total_length;
%w5= length(fitness5') / total_length;
%w6 = length(fitness6') / total_length;
%w7 = length(fitness7') / total_length;

ave_r = [w1 w2 w3 w4] * [r1 r2 r3 r4]';
%ave_r = [w1 w2 w3 w4 w5 w6 ] * [r1 r2 r3 r4 r5 r6 ]';

P = polyfit (all_Energy, all_fitness, 1);
x = -3 :0.5:3;
y = P(1)*x +P(2);

figure
hold on 
plot(Energy1_norm,fitness1_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','cyan', 'LineWidth',line_width)
plot(Energy2_norm,fitness2_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','green', 'LineWidth',line_width)
plot(Energy3_norm,fitness3_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','blue', 'LineWidth',line_width)
plot(Energy4_norm , fitness4_norm, 'o' , 'MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','red', 'LineWidth',line_width)
%plot(Energy5_norm , fitness5_norm, 'o' , 'MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','m', 'LineWidth',line_width)
%plot(Energy6_norm , fitness6_norm, 'o' , 'MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','y', 'LineWidth',line_width)
%plot(Energy7_norm , fitness7_norm, 'o' , 'MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','k', 'LineWidth',line_width)
plot(x,y,'k--','LineWidth',1);
xlabel('Energy (normalized)')
ylabel('Experimental fitness (normalized)')
text(-1.75 ,-1,sprintf('$$r= %.4f $$', ave_r), 'Interpreter','latex', 'FontSize', 16)
legend('Esteban-Riesco2013-1', 'Esteban-Riesco2013-2', 'Esteban-Riesco2013-3', 'Prentoe2019')

 %,'Alzua2022','Serre2013')


if conserve == 1
    title('Conservational-only')
else 
    title('Our Model')
end





