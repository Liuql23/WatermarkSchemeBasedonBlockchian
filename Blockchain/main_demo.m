clear all
close all
clc
% The size of watermark
W_H = 64;W_W = 64;

%���� cost
load cost                                                                                                               
cost_0 = cost{1,5}{1};
[C_H,C_W] = size(cost_0);
[sort_cost,index] = sort(cost_0(:));
aim_coord = index(1:W_H*W_W);
index_mincost = index(1); 


% ���ͼ������С�ĵ�
sign_A = zeros(C_H,C_W);
sign_A(aim_coord)=1;

% ͼ��ֿ�ߴ�
sub_h = 4;
sub_w = 4;
num_aim_block = W_H*W_W/(sub_h*sub_w);
% ����ͼ����������
A = 1:1:C_H*C_W;
A = reshape(A,C_H,C_W);
B_H = C_H-mod (C_H,sub_h);
B_W = C_W-mod (C_W,sub_w);
% ���ɿ��������ֿ���м����
temph = ceil((C_H-B_H)/2);
tempw = ceil((C_W-B_W)/2);
B = A(temph+1:temph+B_H,tempw+1:tempw+B_W);

sign_B = sign_A(temph+1:temph+B_H,tempw+1:tempw+B_W);
% ����sign_B�ֿ�
sign_B_cell=mat2cell(sign_B,sub_h*ones(1,B_H/sub_h),sub_w*ones(1,B_W/sub_w));
sign_B_cell = reshape(sign_B_cell,numel(sign_B_cell),1);
num_block = numel(sign_B_cell);
for i = 1:num_block
    sum_sign_B_cell(i) = sum(sum(sign_B_cell{i}));
end
[sign_B_cell_sort,index_0] = sort(sum_sign_B_cell(:));
sum(sign_B_cell_sort(end-num_aim_block+1:end))

% ����B�ֿ�
B_cell=mat2cell(B,sub_h*ones(1,B_H/sub_h),sub_w*ones(1,B_W/sub_w));
B_cell = reshape(B_cell,numel(B_cell),1);
num_block = numel(B_cell);
init_index = index_0(end-num_aim_block/4-1:end);

% ָ��hash�㷨
alg='MD5'; % hash�㷨{'MD2','MD5','SHA-1','SHA-256','SHA-384','SHA-512'}
sum_time = 0;
T = 0;
while sum_time <= (1/2)*24*60*60
    start = tic;
    % �������һ���ַ���key
    Char_Length = 512; % ����ַ�������
    Limit_Char = [32 127]; %  ����ַ���ȡ��Χ, ��������������ӦASCII��λ��
    key_index = Limit_Char(1) + round((Limit_Char(2)-Limit_Char(1))*rand(1,Char_Length));
    key = char(key_index); % 512λ���ַ���
    
    coord = init_index;
   
    % ������ʼ�飬key + Ini_block_bin
    block_1 = [key,num2str(init_index(end))];
    [H1,index_1]=gen_BlockIndex0(block_1,alg,num_block);
    coord = [coord;index_1];
    
    while numel(coord) < num_aim_block
        % ������һ���� 
        block_next = [key,H1,num2str(index_1)];
        [H1,index_1]=gen_BlockIndex0(block_next,alg,num_block);
        while ismember(index_1,coord)==1
            coord_c = setdiff(1:1:num_block,coord);
            temp_index = ceil(numel(coord_c)/2);
            index_1 = coord_c(temp_index);
        end
        coord = [coord;index_1];
    end
    our_coord = [];
    for coord_i = 1:numel(coord)
        temp = B_cell{coord(coord_i)};
        our_coord = [our_coord;temp(:)];
    end
    time = toc(start);    
    sum_time = sum_time + time;
    count = sum(ismember(our_coord,aim_coord));
%     fprintf('\n count: %d, run time is %.2f seconds.', count, sum_time);
    if count > T
        T = count;
        fprintf('\n count: %d, run time is %.2f seconds.', T, sum_time);
    end
    
end




