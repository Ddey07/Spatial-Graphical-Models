function example()
    %note that, the graph G input must be a connected graph
%=======An Example of the Moral Graph of the Bayesian network Asia=================
    GraphName = 'MoralGraphAsia.txt';%graph G
    %construct a Gaussian distribution N(mu, sigma) that obey  the pairwise Markov
    %property w.r.t the graph G,
    %draw samples from N(mu, sigma) and save the in 'samples.mat' 
    samples = RandomSigmaSamples(GraphName);
%     fname=[GraphName(1:end-4), '_samples',  '.mat'];
%     load(fname, 'samples');
     whos samples

%   Improved IPS procedure    
    [imS, iTer1] = iips(GraphName, samples)
%   The conventional IPS procedure by adjusting the covariance matrix 
    [S, iTer2] = ips1(GraphName, samples)
%   The conventional IPS procedure by adjusting the concentration matrix (inverse of convariance matrix)
    [omega, iTer3] = ips2(GraphName, samples)
%   The Localization Approach proposed by Hara and Takemura (2008)
    [ht, iTer4] = ipsHT(GraphName, samples)
    
    precision(1) = max(max(abs(S-imS)));    
    precision(2) = max(max(abs(omega -ht)));
    precision(3) = max(max(abs(inv(S) -ht)));
        
    maxPrec = max(precision)
%========An Example with the file "graphs/200_0.01_1.txt" input to each procedure==========
    vNum = 200; %the number of vertices in the graph
    prob = 0.01;% the probability of adding edges between two vertices
    i = 1;      %the first graph constructed by¡°constrcutRandomGraph.cpp¡±
    GraphName = ['graphs/', num2str(vNum), '_', num2str(prob), '_', num2str(i), '.txt'];
    samples = RandomSigmaSamples(GraphName);
    whos samples
    tic
    [imS, iTer1] = iips(GraphName, samples);
    whos imS
    iTer1
    toc
    tic
    [S, iTer2] = ips1(GraphName, samples);
    whos S
    iTer2
    toc
    tic
    [omega, iTer3] = ips2(GraphName, samples);
    whos omega
    iTer3
    toc
    tic
    [ht, iTer4] = ipsHT(GraphName, samples);
    whos ht
    iTer4
    toc
    precision(1) = max(max(abs(S-imS)));    
    precision(2) = max(max(abs(omega -ht)));
    precision(3) = max(max(abs(inv(S) -ht)));
    
    maxPrec = max(precision)

        