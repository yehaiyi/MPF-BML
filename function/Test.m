function []= Test(msa_bin,weight_seq,amino_single_combine_array,ind_conserve,J_MPF_BML,con_flag)

markersize = 6;
line_width = 0.5;
% conserved-only model
Single_Mutation_Observed = sum(msa_bin.*weight_seq,1)/sum(weight_seq);  % f_i(a)
H = -log((Single_Mutation_Observed)./(1-Single_Mutation_Observed));

%% Prentoe2019 group4
conserve = con_flag;
[  ~ , test1_seq ] = fastaread('data/Prentoe2019.fasta');
test1_seq = cell2mat(test1_seq');
test1_seq_bin = Binary_Seq(test1_seq,amino_single_combine_array,ind_conserve);

if conserve == 1
    Energy4 = Calculate_Energy(test1_seq_bin, J_MPF_BML,H);
else 
    Energy4 = Calculate_Energy(test1_seq_bin, J_MPF_BML);
end

%Energy4_approx=Energy4;
Energy4_approx = [Energy4(1:2)'  Energy4(3) (Energy4(4)+Energy4(7))/2  Energy4(8) ]' ;
Energy4_norm = normalize(Energy4_approx);


fitness4 = [100, 144.9275, 129.4686, 8.696, 43.478, 133.3333, 10.628, 38.6473, 84.058]' ;
%fitness4_approx =fitness4;
fitness4_approx = [fitness4(1:2)'  fitness4(3) (fitness4(4)+fitness4(7))/2  fitness4(8)]' ; 
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

text(Energy4_norm(1),fitness4_norm(1)+0.1,'S52/JFH1')
text(Energy4_norm(2),fitness4_norm(2)+0.1,'N417D')
text(Energy4_norm(3),fitness4_norm(3)+0.1,'N430D')
text(Energy4_norm(4),fitness4_norm(4)+0.1,'N448D')
%text(Energy4_norm(5),fitness4_norm(5)+0.1,'N476D')
%text(Energy4_norm(6),fitness4_norm(6)+0.1,'N532D')
%text(Energy4_norm(5),fitness4_norm(5)+0.1,'N556D')
%text(Energy4_norm(6),fitness4_norm(6)+0.1,'N623D')
%text(Energy4_norm(9),fitness4_norm(9)+0.1,'N645D')

if conserve == 1
    title('Conservational-only')
else 
    title('Our Model')
end




legend('Prentoe2019')
hold off

end