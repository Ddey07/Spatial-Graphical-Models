function samples = RandomSigmaSamples(GraphName)
%Construct a random Sigma such that N(mu, Sigma) that obey the pairwise Markov
%property w.r.t the graph G,
%Draw samples from N(mu, Sigma) and save the in 'samples.mat' 
%cases is the number of iid samples

    
    cases=5000;%sample size
    %===================input must be a string=======================
    if ~ischar(GraphName)
      error('Input must be a string in the function "input_graph(GraphName)"');
    end
    %===================load a graph============================
    adj = GraphToMatlab(GraphName);
    n = numel(adj);
    adjacencyMatrix=zeros(n);
    for i=1:n
        adjacencyMatrix(i, adj{i})=1; 
    end
   %============costruct a positive matrix=============== 
    global true_cov;
    min_eigen_value=0;
    while min_eigen_value<=0
        while min_eigen_value<=0 
            true_cov=rand(n)*2-1;
            true_cov=(tril(true_cov)+tril(true_cov)').*adjacencyMatrix; 
            for i=1:n
                true_cov(i,i)=sum(abs(true_cov(i,:)))+rand(1);
            end        
            min_eigen_value=min(eig(true_cov));
        end        

        true_cov=inv(true_cov);
        min_eigen_value=min(eig(true_cov));
    end
    clear adjacencyMatrix;
    %=========construct covariance matrix of the graphcial model===========
    %fistly, we change the matrix z and let the diglog
    variance=ceil(rand(1, n)*10);
    for i=1:n
        zz=true_cov(i,i);
        true_cov(i,:)=true_cov(i,:)/sqrt(zz)*sqrt(variance(i));
        true_cov(:,i)=true_cov(:,i)/sqrt(zz)*sqrt(variance(i));
    end

    mu=zeros(1, n);    
    samples=mvnrnd(mu,true_cov,cases);
% %     save('samples.mat', 'samples');
%     fname=[GraphName(1:end-4), '_samples',  '.mat'];
%     save(fname, 'samples');
