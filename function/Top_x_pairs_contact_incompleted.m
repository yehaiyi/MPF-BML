function [toppair_contact ,top_value_real_rank]= Top_x_pairs_contact_incompleted(F_APC,N,x,contact_mat,ind_non_conserve,real_index)

    % 先把real_index转成F_APC里面对应的序号
    F_APC_idx = [];
    for idx = real_index
        idx =idx -383;  % for 4mwf
        idx = find(ind_non_conserve == idx);
        if isempty(idx)    % conserved residue
            continue         
        else
            F_APC_idx = [F_APC_idx idx];
        end
    end

    n = length(F_APC_idx);

    A= reshape(F_APC(F_APC_idx,F_APC_idx),1,n*n);
    top_value = maxk(A,x);

    top_value_real_rank = zeros(1,x);

    for i = 1:x
        top_value_real_rank(i) = length(find(F_APC > top_value(i)))+1;
    end

    toppair_contact = zeros(1,x) -1;


for i = 1:x
       
        [ r_i , c_i ] = find(F_APC== top_value(i));
        I = ind_non_conserve(r_i);
        J = ind_non_conserve(c_i);
        toppair_contact(i) = contact_mat(I,J);
           
end