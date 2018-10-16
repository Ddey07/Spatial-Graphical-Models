% function example()
    clear all
    clc
    %note that, the graph G, i.e., GraphName must be a connected graph
%     GraphName = 'MoralGraphAsia.txt';%graph G
    GraphName = 'graphs\600_0.01_1.txt';%graph G
    %construct a Gaussian distribution N(mu, sigma) that obey  the pairwise Markov
    %property w.r.t the graph G,
    %draw samples from N(mu, sigma) and save the in 'samples.mat' 
    [samples,  true_cov] = RandomSigmaSamples(GraphName);
% %     fname=[GraphName(1:end-4), '_samples',  '.mat'];
% %     load(fname, 'samples');

    nparts = 40;
    threshold = 1e-5;
% %         ===================ips1 and ips2 =================
%         
%         tic
%         s1 = ips1(GraphName, samples);
%         timeS(1) = toc;
% 
%         tic
%         om1 = ips2(GraphName, samples);
%         timeOM(1) = toc;
%       ==================iips and iht=============
        tic
        s2 = iips(GraphName, samples, threshold);%Improved IPS procedure by Xu, Guo and He (2011)
        timeS(2) = toc;

        tic
        om2 = iht(GraphName, samples, threshold);%Improved HT procedure by Xu, Guo and Tang (2012)
        timeOM(2) = toc;
        
 %       ==================ipsp1 and ipsp2=============
        initnparts = nparts;
        sa_input = [2, 1e4, 200, 500, 1e-10, 0];        
        tic
        s3 = ipsp1_rec(GraphName, samples, initnparts, sa_input);
        timeS(3) = toc;
        
        initnparts = nparts;
        tic
        om3 = ipsp2_rec(GraphName, samples, initnparts);
        timeOM(3) = toc;
%       ================== merge ===========================

        partOrNot = -1;
        tic
        sim = chooseSimSA_iipsmp(GraphName, -1);
        time_SAS(4) = toc;     
%         tic
        s4 = iipsmp(GraphName, samples, threshold, sim, partOrNot);
        timeS(4) = toc;
        
        
        partOrNot = -1;
        tic
        sim = chooseSimSA_ihtmp(GraphName, -1);
        time_SAOM(4) = toc;        
%         tic
        om4 = ihtmp(GraphName, samples, threshold, sim, partOrNot);
        timeOM(4) = toc;
%       ======================= merge_partation ===============

        partOrNot = 1;
        tic
        sim = chooseSimSA_iipsmp(GraphName, -1);
        time_SAS(5) = toc
%         tic
        [s, timeSA] = iipsmp(GraphName, samples, threshold, sim, partOrNot);
        timeS(5) = toc
        
        tic
        sim = chooseSimSA_ihtmp(GraphName, -1);
        time_SAOM(5) = toc
%         tic
        om = ihtmp(GraphName, samples, threshold, sim, partOrNot);
        timeOM(5) = toc
        
     %============= precision ===========
%         ps(2) = max(max(abs(s1-s2)));
        ps(3) = max(max(abs(s2-s3)));
        ps(4) = max(max(abs(s2-s4)))
        
%         pom(2) = max(max(abs(om1-om2)));
        pom(3) = max(max(abs(om2-om3)));
        pom(4) = max(max(abs(om2-om4)))
        
%         psom = max(max(abs(s1-inv(om1))))
        

    

        