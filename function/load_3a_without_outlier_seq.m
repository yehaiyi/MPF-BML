function [ ] = load_3a_without_outlier_seq()

load data\NumofPatient_3aE2.mat;
inputfile = 'data/3a_E2_ori.fasta';
[Header_fasta, Sequence_fasta] = fastaread(inputfile);
msa_aa = cell2mat(Sequence_fasta');

% preprocess 
no_patient_idx = find(~patient);
if length(no_patient_idx) >0
    msa_aa(no_patient_idx) = []
end

load data\outliers.mat
msa_aa(outliers,:) =[];
patient(outliers,:) = [];

weight_seq = get_seq_weight(patient);

end