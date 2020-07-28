function [Demos,T,Qi,Qf,Qdi,Qdf,minqPrime1,maxqPrime1,V] = getDemos(e, demos,numDims,numDemos)
Demos = {}; % joint trajectories for exercise e
T = {}; % times for exercise e
Qi = {}; % intial positons
Qf = {}; % final positions
Qdi = {}; % intial velocities
Qdf = {}; % final velocities
Qi{numDims,numDemos} = [];
Qf{numDims,numDemos} = [];
Qdi{numDims,numDemos} = [];
Qdf{numDims,numDemos} = [];
minqPrime1 = 999; % minimum q'_1 in data
maxqPrime1 = -999; % maximum q'_1 in data
for d = 1:length(demos)
    [t, D, Dd, Ddd, Dprime,Ddprime,Dddprime, dt, V] = loadData(e,demos(d));
    if e > 4
        [t, D, Dd, Ddd, Dprime,Ddprime,Dddprime, dt, V] = loadData(e-4,demos(d));
        [t, D, Dd, Ddd, Dprime,Ddprime,Dddprime, dt, V] = loadData(e,demos(d),V);
    end
    demo = D;
    if Dprime(1,1) < minqPrime1
        minqPrime1 =  Dprime(1,1);
    end
    if Dprime(1,end) > maxqPrime1
        maxqPrime1 =  Dprime(1,end);
    end
    Demos{d} = demo;
    for i = 1:numDims
        Qi{i,d} = D(i,1);
        Qdi{i,d} = Dd(i,1);
        Qf{i, d} = D(i,end);
        Qdf{i,d} = Dd(i,end);
    end
    T{d} = t;
    %     figure(100)
    %     plot(demo(1,:),demo(2:end,:));hold on
    %     figure(101)
    %     plot(Dprime(1,:),Dprime(2:end,:));hold on
    %     figure(102)
    %     plot(Dprime(1,:),.5*Ddprime(1,:).^2);hold on
        figure(103)
        plot(t,D');hold on
        figure(104)
        plot(t,Ddprime');hold on
end
end