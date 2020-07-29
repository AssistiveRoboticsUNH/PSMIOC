%% train models
clear all
close all
addpath GMM-GMR-DMP\
for e = 1:4
    X = [];
    close all
    [Demos,T,Qi,Qf,Qdi,Qdf,minqPrime1,maxqPrime1,V] = getDemos(e, 1:3,4,3);
    for i =1:length(Demos)
        D = Demos{i};
        t = T{i};
        t = t./t(end); % scale time from 0 to 1
        tmp = [t;D];
        X = [X tmp];
    end
    nbStates = 4; % number of components used in GMM.
    Data = X;
    nbVar = size(Data,1);
    % Training of GMM by EM algorithm, initialized by k-means clustering.
    [Priors, Mu, Sigma] = EM_init_kmeans(Data, nbStates);
    [Priors, Mu, Sigma] = EM(Data, Priors, Mu, Sigma);
    expData(1,:) = linspace(min(Data(1,:)), max(Data(1,:)), 200); % expected distribution
    [expData(2:nbVar,:), expSigma] = GMR(Priors, Mu, Sigma,  expData(1,:), [1], [2:nbVar]);
    save(['models/expData_' num2str(e)],'expData')
    save(['models/expSigma_' num2str(e)],'expSigma')
end
%% Experiment 2: patient evaluation
close all
plot(expData(1,:),expData(2:5,:)','LineWidth',2);hold on
plot(X(1,:),X(2:end,:)','Color','k') % plot data and GMM result
expdt = (expData(1,end) - expData(1,1))/length(expData(1,:));
table = zeros(4,8,3);
for e = 1:4
    load(['models/expData_' num2str(e)])
    load(['models/expSigma_' num2str(e)],'expSigma')
    for j = [0 4]
        errorTot = 0;
        for i = 4:6
            [Demos,T,Qi,Qf,Qdi,Qdf,minqPrime1,maxqPrime1,V] = getDemos(e+j, i,4,1);
            D = Demos{1};
            t = T{1};
            t =  t./t(end); % scale time from 0 to 1
            delta = 0*D;
            H = D*0;
            H2 = zeros(1,size(expData,2));
            for i2 = 1:length(t)
                ind = floor(t(i2)/expdt)+1;
                if ind > size(expData,2)
                    ind = size(expData,2);
                end
                delta(:,i2) = (D(:,i2) - expData(2:5,ind));
                H2(ind) = H2(ind)+delta(:,i2)'*inv(expSigma(:,:,ind))*delta(:,i2); % weight offset by covariance matrix
            end
            error = sum(H2);
            table(e,e+j,i-3) = error;
        end
    end
end
figure
plotTable(log(table./1E3))
sum(mean(table,3),1)