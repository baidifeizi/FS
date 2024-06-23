function [outA,outB,outObj,outNumIter] = RSSFCA(X,gama,maxIter,cluster_n)
% A (cxn) is the membership matrix.
% B (mxc) is the cluster center matrix.
c=cluster_n; %设置聚类中心数量 c。
n=size(X,2); %获取数据矩阵 X 的列数 n
thresh = 10^-5; %设置收敛阈值 thresh。
obj = zeros(maxIter,1); %初始化目标函数值数组 obj。
[dim,data_n]=size(X); %获取数据矩阵 X 的维度 dim 和数据点数量 
%设置 FCM 算法的参数。
options = [2;	% exponent for the partition matrix U 分区矩阵 U 的指数
    100;	% max. number of iteration 最大迭代次数
    1e-5;	% min. amount of improvement 最小改进量
    0];	% info display during iteration  在迭代过程中显示信息

[center, A] = fcm(X', c, options); % Initialize parameters ：使用 FCM 算法对数据进行初始聚类，返回聚类中心 center 和隶属度矩阵 A。
B=center'; %将聚类中心矩阵转置为 B。
GMM_WIDTH=1; %设置 GMM（高斯混合模型）的宽度。
%% compute variance 计算方差
for i = 1:c
    diffs = X'-(ones(data_n, 1) * B(:,i)'); %计算数据与聚类中心 i 的差异。
    diffs = diffs.*(sqrt(A(i,:))'*ones(1, dim));  % 对差异加权，基于隶属度矩阵 A。
    covars(:,:,i)=(diffs'*diffs)/sum(A(i,:)); %计算协方差矩阵
    
    %处理协方差矩阵的秩问题
    try   
        if rank(covars(:,:,i)) < dim
            covars(:,:,i) = covars(:,:,i) + GMM_WIDTH.*eye(dim);
        end
    catch
        covars(:,:,i) = covars(:,:,i) + GMM_WIDTH.*eye(dim);
    end
end
%计算距离和更新
%循环遍历每个聚类中心。
for i = 1:c
    diffs = X'-(ones(data_n, 1) * B(:,i)'); %计算数据与聚类中心 i 的差异。
    mahalanobis2 = sum((diffs*inv(covars(:,:,i))) .* diffs,2);   %计算马氏距离平方。
    D(:,i) =1/2.*(mahalanobis2+log(det(covars(:,:,i)))+dim.*log(2*pi));  %计算高斯密度对数。
end


MIN_COVAR =1e-3;  %设置协方差矩阵的最小值。
init_covars=covars; %存储初始协方差矩阵。
for t = 1:maxIter %主循环，迭代进行聚类。
    %t
    %-------- update A when fixed B 当固定 B 时更新隶属度矩阵  A--------
    A = updateA_capped_robust(X,B,gama,D); %A denotes degree of membership--c*n,
    
    %-------- update B when fixed A and lamda 当固定 A 和 lamda 时更新聚类中心 B--------%
    B = updateB_capped_robust(X,A,B); %B denotes centers
    
    %------- calculate the covariance 计算协方差 --------%
    %再次计算协方差矩阵
    for i = 1:c
        diffs = X'-(ones(data_n, 1) * B(:,i)');
        diffs = diffs.*(sqrt(A(i,:))'*ones(1, dim));
        covars(:,:,i)=(diffs'*diffs)/(sum(A(i,:))+eps); %计算协方差矩阵。
    end
    
    P=sum(covars,3);%计算协方差矩阵总和
    Nan=sum(isnan(P(:)));%检查协方差矩阵是否包含 NaN。
    if isreal(covars)==0 || Nan,break; end,%如果协方差矩阵不为实数或包含 NaN，则跳出循环。
    for j = 1:c %遍历协方差矩阵
        if min(svd(covars(:,:,j))) < MIN_COVAR  %检查协方差矩阵的最小值。
            covars(:,:,j) = init_covars(:,:,j); %恢复初始协方差矩阵。
        end
    end

    %计算目标函数和距离
    Pi=sum(A,2)./data_n; %计算每个聚类中心的权重。
    for i = 1:c %计算距离
        diffs = X'-(ones(data_n, 1) * B(:,i)'); %计算数据与聚类中心 i 的差异。
        mahalanobis2 = sum((diffs*inv(covars(:,:,i))) .* diffs,2); %马氏距离平方。
        Temp(:,i) =1/2.*(mahalanobis2+log(det(covars(:,:,i)))+dim.*log(2*pi));  %计算高斯密度对数。
    end
    D=Temp; %更新距离矩阵。
    
    
    for ii=1:n %处理距离矩阵。
        AI=D(ii,:); %获取第 ii 行的距离。
        if sum(AI<0)>1 ||sum(AI<0)==1 %检查是否有负值。
            TempT(ii,:)=AI-min(AI)+eps; %调整距离矩阵。
        else
            TempT(ii,:)=AI;
        end
    end

    obj(t)=sum(sum(TempT'.*A + A.^2*gama)); %计算目标函数。
    
    CC{t}=B; %存储当前聚类中心。您也可以使用目标函数作为收敛条件 You can also use objective funtion as the convergence condition  
    if(t > 1)
        diff = max(max(CC{t-1}-CC{t}));
        if(diff < thresh) %判断是否收敛。
        end
    end
end

% figure,mesh(reshape(E,rows,cols))
outNumIter = t; %实际迭代次数
outObj = obj(1:t);%目标函数值数组。
outA = A; %隶属度矩阵
outB = B;%聚类中心矩阵