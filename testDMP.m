clear all
addpath GMM-GMR-v2.0\
b = 10;
nbStates = 20;
nbVar = 2;
%% learn trajectories
for d = 1:4
    for e = 3
        X = [];
        close all
        for i = 1:3
            [t,D,Dd, Ddd,Dprime,Ddprime,Dddprime, dt,V]  = loadData(e,i);
            D = D(d,:);
            [x,~,s,f] = DMP(t, D,b, D(end),0,0);
            plot(s,f); hold on
            tmp = [s; f];
            X = [X tmp];
        end
        expDMP = [s;f];
        [Priors, Mu, Sigma] = EM_init_kmeans(X, nbStates);
        [Priors, Mu, Sigma] = EM(X, Priors, Mu, Sigma);
        expData(1,:) = flip(linspace(min(X(1,:)), max(X(1,:)), 500));
        [expData(2:nbVar,:), expSigma] = GMR(Priors, Mu, Sigma,  expData(1,:), [1], [2:nbVar]);
        plot(expData(1,:),expData(2,:),'Color','k')
        expDMP = expData(1:2,:);
        save(['expDMP_' num2str(e) '_' num2str(d)],'expDMP')
    end
end
%% test
for e = 3
    figure(e)
        for d = 1:4
            load(['expDMP_' num2str(e) '_' num2str(d)],'expDMP')
            for q = 1:3
                [t,D,Dd, Ddd,Dprime,Ddprime,Dddprime, dt,V]  = loadData(e,q);
                goal = D(:,end);
                subplot(1,4,d)
                plot(t, D(d,:),'Color','k','LineWidth',1);hold on
                [t,D,Dd, Ddd,Dprime,Ddprime,Dddprime, dt,V]  = loadData(e,q);
                [T,y] = DMP(t, D(d,:), b, goal(d), 0, 0,expDMP);
                plot(T, y, 'Color','b','LineWidth',2);
                axis tight
                xlim([0 1])
            end
        end
end
figure(e)
set(gcf,'Color','w')