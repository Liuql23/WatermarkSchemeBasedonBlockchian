function  [H,index_0,index_1]=gen_BlockIndex0(block,alg,num_block)
sup = ceil(log2(num_block));
%�����ʼ���hashֵ
H=hash(block,alg);
%��Hashֵת��Ϊ�����Ʒ�����ȡ��sup�������Ʒ�
H0_bin=dec2bin(H);
bin1 = H0_bin(end-sup+1:end);
bin2 = H0_bin(end-sup-sup+1:end-sup);
% ��binת��Ϊʵ��
index_0_temp = rem(bin2dec(bin1),num_block);
index_1_temp = rem(bin2dec(bin2),num_block);
if index_0_temp == 0
    index_0 = num_block;
else
    index_0 = index_0_temp;
end

if index_1_temp == 0
    index_1 = num_block;
else
    index_1 = index_1_temp;
end

if index_1 == index_0
    if index_0 == num_block
        index_1 = index_1 -1;
    else
        index_1 = index_1+1;
    end
end

end