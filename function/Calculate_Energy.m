function E = Calculate_Energy (seq_bin, J_MPF_BML)

E = diag( seq_bin * J_MPF_BML * seq_bin');
