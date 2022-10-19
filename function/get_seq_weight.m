function weight_seq = get_seq_weight(patient)

idx=1;
weight_seq = zeros([length(patient),1]);
for id=patient'
    weight = 1/ (sum(patient(:)==id));
    weight_seq(idx)=weight;
    idx=idx+1;
end