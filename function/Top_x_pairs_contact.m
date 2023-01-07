function toppair_contact = Top_x_pairs_contact(F_APC,N,x,contact_mat)

    A= reshape(F_APC,1,N*N);
    top_value = maxk(A,x);

    toppair_contact = zeros(1,x);

    for i = 1:x
        [ r_i , c_i ] = find(F_APC== top_value(i));
        toppair_contact(i) = contact_mat(r_i,c_i);
    end

end