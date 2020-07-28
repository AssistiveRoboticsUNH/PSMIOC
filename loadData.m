function [t,D,Dd, Ddd,DPrime,DdPrime,DddPrime, dt,V] = loadData(e,r,Vin)
% e 1-4 corresponds to correct exercises and 5-8 corresponds to incorrect exercises
offset = 0;
if (e==4 || e==8)
    offset = 2;
end

if e<=4
    load(['data/TPR Therapist Weights000' num2str(e+offset) 'Q.mat'])
else
    load(['data/TPR Therapist Weights Bad Compensation000' num2str(e-4+offset) 'Q.mat'])
end
dt = 1/300;
time = size(D{r},2)*dt;
dataLength = 500;
dt = time/dataLength;
numDims = 4;
newData = zeros(numDims,dataLength,length(D));
tf = zeros(1,length(D));
for i = 1:length(D)
    datai = D{i};
    traj = D{i};
    traj = traj(1:numDims ,:);
    newDatai = zeros(numDims,dataLength);
    t = (1:size(datai,2))*dt;
    tf(i) = t(end);
    newT = linspace(1,length(t),dataLength)*dt;
    for j =1:numDims
        newDatai(j,:) = interp1(t,traj(j,:), newT, 'spline', 'extrap');
    end
    newData(:,:,i) = newDatai;
end
D = newData(:,:,r);
D = D(:,10:floor(1/2*dataLength));
D1 = newData(:,:,1);
D1 = D1(:,10:floor(1/2*dataLength));
% allData = newData(:, 1:dataLength/2,:);
if exist('Vin') ~= 1
    V = getV(D1-D1(:,end));%; % this is the eigenvector matrix for a signle demonstration
else
    V = Vin;
end
DPrime = inv(V)*D;
t = (1:size(D,2))*dt;
n = length(t);
h = .06;
for i = 1:size(DPrime,1)
    ri = ksr(t,DPrime(i,:),h,n);
%     DPrime(i,:) = ri.f;% smooth(DPrime(i,:),.2); % smooth data
end
D = V*DPrime;
Dd = 0*D;
Ddd = 0*D;
DdPrime = 0*DPrime;
DddPrime = 0*DPrime;
for i =1:size(D,1)
    Dd(i,:) = d_dt(D(i,:),dt);
    DdPrime(i,:) = d_dt(DPrime(i,:),dt);
    Ddd(i,:) = d_dt(Dd(i,:),dt);
    DddPrime(i,:) = d_dt(DdPrime(i,:),dt);
end
if DPrime(1,end)-DPrime(1,1)<0
    DPrime(1,:) = -DPrime(1,:);
end
for i = 1:size(D,1)
    DdPrime(i,:) = d_dt(DPrime(i,:),dt);
    DddPrime(i,:) = d_dt(DdPrime(i,:),dt);
end
offset = floor(size(DdPrime,2)/2);
% crop motion once velocity falls below a threshold
ind2 = offset+find( (.5*DdPrime(1,offset:end).^2) < 0.1,1)-1;
ind1 = find(.5*DdPrime(1,1:end).^2 > 0.1,1);
if isempty(ind1)
    ind1 = 2;
end
if isempty(ind2)
    ind2 = size(DPrime,2)-1;
end
tmp = ind1:ind2;
DPrime = DPrime(:,tmp);
DdPrime = DdPrime(:,tmp);
DddPrime = DddPrime(:,tmp);
D = D(:,tmp);
Dd = Dd(:,tmp);
Ddd = Ddd(:,tmp);
t = t(:,tmp);
t = t-t(1);

n = length(t);
h = .05;
ri = ksr(t,DPrime(1,:),h,n);
DPrime(1,:) = ri.f;% smooth(DPrime(i,:),.2); % smooth data
h = .2;
for i = 2:size(DPrime,1)
    ri = ksr(DPrime(1,:),DPrime(i,:),h,n);
    DPrime(i,:) = ri.f;% smooth(DPrime(i,:),.2); % smooth data
end
for i = 1:size(D,1)
    DdPrime(i,:) = d_dt(DPrime(i,:),dt);
    DddPrime(i,:) = d_dt(DdPrime(i,:),dt);
end
D = V*DPrime;
Dd = V*DdPrime;
Ddd = V*DddPrime;
offset = floor(size(DdPrime,2)/2);
% crop motion once velocity falls below a threshold
ind2 = offset+find( (.5*DdPrime(1,offset:end).^2) < 0.1,1)-1;
ind1 = find(.5*DdPrime(1,1:end).^2 > 0.1,1);
if isempty(ind1)
    ind1 = 2;
end
if isempty(ind2)
    ind2 = size(DPrime,2)-1;
end
tmp = ind1:ind2;
DPrime = DPrime(:,tmp);
DdPrime = DdPrime(:,tmp);
DddPrime = DddPrime(:,tmp);
D = D(:,tmp);
Dd = Dd(:,tmp);
Ddd = Ddd(:,tmp);
t = t(:,tmp);
t = t-t(1);
end