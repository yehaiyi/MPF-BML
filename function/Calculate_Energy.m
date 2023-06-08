function E = Calculate_Energy (seq_bin, J_MPF_BML,H)

switch nargin 
    case 2
        E = diag( seq_bin * J_MPF_BML * seq_bin');
    case 3 
        E = seq_bin * H' ;
end

end

