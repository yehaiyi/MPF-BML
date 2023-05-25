% read S52 strain with residues from 384 to 752
[  ~ , seq ] = fastaread('data/Yu2013.fasta');
msa = cell2mat(seq)
%% 
[  ~ , seq ] = fastaread('data/preprocess/GT3a.D17763_384-752.fasta');

res=583

if res>479 && res<579
    res=res+1;
else if res>579
    res=res+6;
end
end

seq(res-383)
%% 
[  ~ , seq ] = fastaread('data/preprocess/GT3a.D17763_384-752.fasta');



%%  比较在E2上突变的数量
clear
[  ~ , Urban ] = fastaread('data/Prentoe2019.fasta');
[  ~ , seq ] = fastaread('data/preprocess/S52.JFH1_EU204645.fasta');
S52=seq(384:752);

for i = 1:length(Urban)
    seqi =  cell2mat(Urban(i));
    num_mutant = sum(~(S52==seqi))
    find(~(S52==seqi));
    
end

%% 
clear
[  ~ , seq ] = fastaread('data/preprocess/S52.JFH1_EU204645.fasta');

for res =[417 430 448 476 532 556 623 645 ]
    res
    if res>479 && res<579
        res=res+1;
    else if res>579
        res=res+6;
    end
    end

    if seq(res)=='N'
        seq(res)='D';
    else 
        error("!")
    end

    seq(384:752)
    [  ~ , seq ] = fastaread('data/preprocess/S52.JFH1_EU204645.fasta');
end 


%% Alignment

[  ~ , seq ] = fastaread('data/preprocess/GT3a.D17763_compare.fasta');
Aseq = multialign(seq)
%% 特定位点上的出现过的aa及频率计算（要算上patient_weight的）
%这个也得是没剔除conserved res的msa_aa ( i.e. 使用msa_aa_369)
% 导入369的msa_aa
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
load data\weight_seq.mat
msa_aa_369(outliers,:) =[];
patient(outliers,:) = [];
%
res=408

if res>479 && res<579
    res=res+1;
else if res>579
    res=res+6;
end
end

aas = msa_aa_369(:,res-383);
aa_unique = unique(aas);
num_patient = sum(weight_seq);
freq_list = zeros(length(aa_unique),1);
for i =1: length(aa_unique)
    aa = aa_unique(i);
    freq_list(i) =  sum(weight_seq .* (aas==aa) )  / num_patient;
end

table(aa_unique,freq_list)
%sortrows(table,'freq_list', 'descend')

index_after = find(ind_non_conserve== res-383)
c=amino_single_combine_array(index_after);
c{1,1}

%%
%特定位点上mutant combine之后的aa list
% 使用前要导入模型那里
res=417;

if res>479 && res<579
    res=res+1;
else if res>579
    res=res+6;
end
end

index_after = find(ind_non_conserve== res-383);
c=amino_single_combine_array(index_after);
c{1,1}

%%  找出conserved residues
% 这里的msa_aa需要是还没有剔除conserved res的，所以num_res应该要=369
num_res= size(msa_aa,2);
num_mutant_list = zeros(1,num_res);

for i = 1:num_res
    aas = msa_aa(:,i);
    num_mutant_list(i) = length(unique(aas));
end

conserved_idx = find(num_mutant_list==1);
length(conserved_idx)

%% 原始residue index与去掉conserved residue之后的index转换
res=522;

if res>479 && res<579
    res=res+1;
else if res>579
    res=res+6;
end
end

index_before=res-383;
i=0;
while(index_before>(conserved_idx(i+1)) && i <=length(conserved_idx) )
    i=i+1;
end

index_after=index_before-i;
amino_single_combine_array(index_after)

%% 单突变对应的binary_seq
[  ~ , S52 ] = fastaread('data/preprocess/S52.JFH1_EU204645.fasta');
S52_seq = S52(384:752);
mut_seq=mutate(S52_seq,645,'D');

S52_bin_seq= Binary_Seq(S52_seq,amino_single_combine_array,ind_conserve);
mut_bin_seq=Binary_Seq(mut_seq,amino_single_combine_array,ind_conserve);

  
com_seq = S52_bin_seq - mut_bin_seq;
idx_change = find(com_seq)   
value_change = com_seq(idx_change)
mut_bin_ind = find(com_seq== -1);

%% 某个突变对应的Energy变化效果
column_idx = 546

h_effect =J_MPF_BML(column_idx,column_idx)
J=J_MPF_BML(:,column_idx);
background = bin_seq;
background(column_idx)=0;
J_effect = J .* background';
J_effect=J_effect(find(J_effect))
non_zero_J = length(J_effect)
Energy_change = h_effect +2* sum(J_effect)
h_conserved =H(column_idx)


%% bin index所对应的res
bin_idx = 25
sum_aa = 0;
flag=0;
res_idx_aft=1;

while(flag==0)
    num_aa = length(amino_single_combine_array{res_idx_aft,1}) -1;
    sum_aa = sum_aa + num_aa;
    if bin_idx> sum_aa
        res_idx_aft=res_idx_aft+1;
    else 
        flag=1;
    end
end

res_idx_aft
rank = sum_aa - bin_idx +2 
amino_single_combine_array{res_idx_aft,1}
%res_idx_bf = ind_non_conserve(res_idx_aft) + 383 % S52的numbering

%% 将msa_aa中的gap替换成真实的aa（按概率抽样）

s = RandStream('mrg32k3a');


gap_col_idx = 1;

aas = msa_aa(:,gap_col_idx);
aa_unique = unique(aas);
num_patient = sum(weight_seq);
freq_list = zeros(length(aa_unique),1);
for i =1: length(aa_unique)
    aa = aa_unique(i);
    if aa=='-'
        freq_list(i)=0;
    else
        freq_list(i) =  sum(weight_seq .* (aas==aa) )  / num_patient;
    end
end

freq_list = freq_list / sum(freq_list);
alphabet = 1:length(aa_unique);
prob = freq_list;
out = randsrc(1,1,[alphabet;prob'],s);
new_aa = aa_unique(out)
%% 
new_msa = GapSubstitution(msa_aa,weight_seq);
save data\msa_without_gaps.mat new_msa
%% 
ind_conserve_H77 = [];

for index_before =  ind_conserve +383

    if index_before>480 && index_before<580
        index_before = index_before - 1;
    else if index_before > 585
            index_before =index_before -6;
    end
    end
    ind_conserve_H77 = [ind_conserve_H77 index_before];

end

%% 4mwf 的index转换
real_index =[];

for index = location_index

    if index>=1 && index <= 48
        index = index + 411;
        real_index =[real_index index];

    elseif index >=49 && index <= 212
        index = index + 433;
        real_index =[real_index index];
    end
end


        
        





















