% model - metabolic model in standard format
% exp_mat - expression matrix in which row j and column i has [e^j_i]/[e^ref_i]
% met_frw_mat - saturation values matrix in which row j and column i has a+^j_i
% met_bck_mat - saturation values matrix in which row j and column i has a-^j_i
% rxns_idx - A vector with  reactions indices corresponding to the core set
% used in the problem
% ko_rxns - A matrix where for each row (knockout j), the indices of the
% knocked-out reactions are specify
% ex_flux - exchange fluxes matrix where each row represent a different exchange reaction
% and each column represent a differen knockout
% ex_rxns - A vector with the corresponding exchange reaction indices 
% GR - A vector of growth rates values in each one of the knockouts

function [V_v] = IOMA(model, exp_mat, met_frw_mat,met_bck_mat,rxns_idx,ko_rxns,ex_flux,ex_rxns,GR)

[met_num, rxn_num] = size(model.S);
[ko_exp_num] = size(exp_mat, 1);
biomass_idx = find(model.c==1);

reactions_with_v_max = find(sum(met_frw_mat~=0));
v_max_num = length(reactions_with_v_max);

ko_num = ko_exp_num*ones(1,v_max_num);
MAX_FLUX = 100;
MIN_FLUX = -100;

%QP variables (v_max_frw, v_max_bck, v1, e1,y1 ..vk, ek,yk)
variables_number = 2*v_max_num  +  ko_exp_num*(rxn_num + 2*v_max_num);
ko_vars =  (rxn_num + 2*v_max_num);
m = [];
temp_size1 = 0;
temp_size2 = 0;
for i=1:ko_exp_num
    %Constructing Mass Balance constraints matrix
    m_zeros =  sparse(met_num,1*temp_size1);
    m_b_constraints = sparse(met_num, 2*v_max_num);
    m_b_constraints = [m_b_constraints,m_zeros, model.S, sparse(met_num, 2*v_max_num), sparse(met_num, (ko_exp_num-i)*ko_vars)];
    
    %Constructing Michaelis Menten constraints matrix
    exp_e = exp_mat(i, :);
    exp_ce = exp_mat(i, :) .* met_frw_mat(i, :);
    exp_de = exp_mat(i, :) .* met_bck_mat(i, :);
    one_vec = ones(v_max_num,1);
    one_vec(nan_idx) = 0;
    
    m2f = -sparse([ 1:v_max_num],[ 1:v_max_num], exp_ce); %v_max+
    m2b = sparse([ 1:v_max_num],[ 1:v_max_num], exp_de); %v_max-
    m_zeros2 =  sparse(v_max_num,1*temp_size1);
    m3 = [sparse([1:v_max_num],rxns_idx, one_vec), sparse(v_max_num,(rxn_num-max(rxns_idx)))];%v
    m4 = -sparse([1:v_max_num],[1:v_max_num], exp_e);%e
    m5 = sparse(v_max_num, v_max_num);%y
    MM_constraints = [m2f,m2b, m_zeros2, m3, m4, m5, sparse(v_max_num, (ko_exp_num-i)*ko_vars)];
    temp_size1 = temp_size1 + (rxn_num + 2*v_max_num);
    
    %Constructing the variance constraints matrix
    e1_vec = ones(v_max_num,1).*((ko_num'-ones(v_max_num,1))./(ko_num'));%(k-1)/k
    e2_vec = (ones(v_max_num,1)./(ko_num'));%(1/k)
    y_vec = sqrt(ko_exp_num)*ones(v_max_num,1);
    m_zeros3 = sparse(v_max_num, v_max_num);%zeros
    m_vmax = sparse(v_max_num, 2*v_max_num);%v_max+&-
    m_s = sparse(v_max_num, rxn_num);%Stoichiometric matrix
    m_e1 = -sparse([1:v_max_num],[1:v_max_num], e1_vec);%e1
    m_y = sparse([1:v_max_num],[1:v_max_num], y_vec);%y
    m_e2 = sparse([1:v_max_num],[1:v_max_num], e2_vec);%e2
    m_temp1 = [m_s, m_e1, m_y];
    m_temp2 = [m_s, m_e2, m_zeros3];
    QP_constrains = [m_vmax];
    for j=1:temp_size2
        QP_constrains = [QP_constrains, m_temp2];
    end
    QP_constrains = [QP_constrains, m_temp1];
    j = ko_exp_num-i;
    while j>=1
        QP_constrains = [QP_constrains, m_temp2];
        j = j - 1;
    end
    temp_size2 = temp_size2 + 1;
    
    %combining all three matrices
    m = [m;m_b_constraints; MM_constraints; QP_constrains];
end

%Setting bounds and indices
count_y = rxn_num + 3*v_max_num;%y variables
count_v = 2*v_max_num;%v variables
count_e = rxn_num + 2*v_max_num;%epsilon variables
count_b = 2*v_max_num;%biomass variables

ub = zeros(variables_number,1);
lb = zeros(variables_number,1);
%v_max bounds
lb(1:v_max_num) = 0;
ub(1:v_max_num) = MAX_FLUX;
lb(v_max_num+1:2*v_max_num) = MIN_FLUX;
ub(v_max_num+1:2*v_max_num) = 0;

%v bounds
model.lb(find(model.lb == -1000)) = MIN_FLUX;
model.ub(find(model.ub == 1000)) = MAX_FLUX;

y_idx = [];
v_idx = [];
e_idx = [];
b_idx = [];

for i=1:ko_exp_num
    y_idx = [y_idx, [count_y+1:count_y+v_max_num]];
    e_idx = [e_idx, [count_e+1:count_e+v_max_num]];
    v_idx = [v_idx, [count_v+1:count_v+rxn_num]];
    b_idx = [b_idx, [count_b+biomass_idx]];
    
    %v bounds
    ub(count_v+1:count_v+rxn_num) = model.ub;
    lb(count_v+1:count_v+rxn_num) = model.lb;
    %exchange reaction bounds
    ub(ex_rxns'+count_v.*ones(length(ex_rxns),1)) = ex_flux(:,i);
    lb(ex_rxns'+count_v.*ones(length(ex_rxns),1)) = ex_flux(:,i);
    
    %set KO bounds
    ko = ko_rxns(i,:);
    ko = ko(ko~=0);
    ub(ko'+count_v.*ones(length(ko),1)) = 0;
    lb(ko'+count_v.*ones(length(ko),1)) = 0;
    
    %update counters
    count_y = count_y + rxn_num + 2*v_max_num;
    count_v = count_v + rxn_num + 2*v_max_num;
    count_e = count_e + rxn_num + 2*v_max_num;
    count_b = count_b + rxn_num + 2*v_max_num;
end

%epsilon, y and GR bounds
lb(e_idx) =-Inf;
ub(e_idx) = Inf;
lb(y_idx) = -Inf;
ub(y_idx) = Inf;
lb(b_idx) = GR;
ub(b_idx) = GR;

%prepare problem for tomlab
Name = model.description;
A = m;
F = sparse(y_idx',y_idx', 2*ones(ko_exp_num*v_max_num,1));
x_L = lb;
x_U = ub;
b_U = zeros(size(m,1),1);
b_L = zeros(size(m,1),1);
c = zeros(variables_number,1);
[m,n] = size(A);
x_0 = zeros(n,1);
x_min = x_L; x_max = x_U;
fprintf('qp problem. Variables %d. Knapsacks %d\n',n,m);
Prob = qpAssign(F, c, A, b_L, b_U, x_L, x_U, x_0, Name);
Prob.MIP.cpxControl.QPMETHOD = 1;

model.tomlab_result = tomRun('cplex', Prob, 1);
model.result_vector = model.tomlab_result.x_k;
model.result_opt = model.tomlab_result.f_k;
model.result_status = model.tomlab_result.ExitFlag;
model.result_status_text = model.tomlab_result.ExitText;
v_result = model.result_vector(v_idx);
V_v = reshape(v_result,rxn_num,ko_exp_num);
fprintf('\n*** RunTomlabQP ***\nOpt val: %d\nExit flag: %d\nExit text: %s\n',model.result_opt,model.result_status,model.result_status_text);



