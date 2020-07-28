function [X,Y,sList,f] = DMP(t,xn,b,goal,offset,offsetv, expDMP)
nbData = length(xn); %Length of each trajectory
nbVar = 1; %Number of variables (Trajectory in a plane)
nbStates = b; %Number of states (or primitives)
nbRepros = 1;
kP = 50; %Initial stiffness gain
kV = 2*sqrt(kP);
% kP = 100000; %Initial stiffness gain
% kV = 20000; %Damping gain
dt = (t(2)-t(1)); %Time step
%Decay factors
alpha = 1;
%Centers equally distributed
Mu_d = linspace(nbData,1,nbStates);
Sigma_d = 100;
%Estimate Mu_s and Sigma_s to match Mu_d and Sigma_d
Mu_s(1,:) = exp(-alpha*Mu_d*dt);
for i=1:nbStates
    std_s = Mu_s(1,i) - exp(-alpha*(Mu_d(i)+Sigma_d^.5)*dt);
    Sigma_s(:,:,i) = std_s^2;
end
alpha = 1;
%% format
Data = zeros(2,length(xn));
Data(1,:) = xn';
% Data(2,:) = yn';
posId = [1:nbVar]; velId=[nbVar+1:2*nbVar]; accId = [2*nbVar+1:3*nbVar];
%Define end-point as target
xT = Data(posId,end);
%Estimate of derivatives with polynomial fitting
Data(velId,:) = computeDerivative(Data(posId,:), dt); %Velocity
Data(accId,:) = computeDerivative(Data(velId,:), dt); %Acceleration

%% Learning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WLS learning: least norm solution to find Mu_X (Y=Mu_x*H')
s = 1; % Decay termY
S = [];
for n=1:nbData
    s = s + (-alpha*s)*dt;
    sList(n) = s;
    for i=1:nbStates
        h(i) = gaussPDF(s,Mu_s(:,i),Sigma_s(:,:,i));
        
    end
    %Compute weights
    H(n,:) = h./sum(h);
end
Y = (Data(accId,:)-kP*(Data(posId,end) - Data(posId,:))+kV*(Data(velId,:)))./(sList ); % ((Data(posId,end) - Data(posId,1)))
Mu_F = [inv(H'*H)*H'*Y']';
% plot((H'.*Mu_F')')
%% Reproductions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in=[1]; out=[2:3];
xT = goal;
bool = 0;
if exist('expDMP')==1
    bool = 1;
end
for nb=1:nbRepros
    if nb==1
        rWLS(nb).currPos = Data(posId,1)+offset; %Initial position
        rWLS(nb).currVel = zeros(nbVar,1)+Data(velId,1)+offsetv; %Initial velocity
        rWLS(nb).currAcc = zeros(nbVar,1)+Data(accId,1)*0; %Initial acceleration
    end
    s = 1; %Decay term
    n2 = 1;
    f = zeros(1,length(nbData));
    n = 0;
    n3 = 0;
    %   for n = 1:nbData
    % while abs(rWLS(nb).currPos-goal)>0.001 || n < nbData
    while n3 < nbData
        n = n+1;
        if n==nbData+1
            n = n-1;
            n3 = n3+1;
        end
        %Log data
        rWLS(nb).Data(:,n+n3) = [rWLS(nb).currPos; rWLS(nb).currVel; rWLS(nb).currAcc];
        s = s + (-alpha*s)*dt;
        for i=1:nbStates
            rWLS(nb).H(i,n) = gaussPDF(s,Mu_s(:,i),Sigma_s(:,:,i));
        end
        rWLS(nb).H(:,n) = rWLS(nb).H(:,n) /sum(rWLS(nb).H(:,n));
        rWLS(nb).currAcc = kP*(goal- rWLS(nb).currPos)-kV*rWLS(nb).currVel+Mu_F*rWLS(nb).H(:,n)*(s); % *(goal - Data(posId,1))
        if n3 == 0
            f(n) = Mu_F*rWLS(nb).H(:,n);
            sList(n) = s;
        end
        if bool==1
            while s < expDMP(1,n2) && n2< size(expDMP,2)
                n2 = n2+1;
            end
            if n2 >= size(expDMP,2)
                expDMP(2,n2) = 0;
            end
            rWLS(nb).currAcc = kP*(goal- rWLS(nb).currPos)-kV*rWLS(nb).currVel+expDMP(2,n2)*(s); % *(goal - Data(posId,1))
        end
        % rWLS(nb).currAcc = (xT-rWLS(nb).currPos)*kP - rWLS(nb).currVel*kV+Mu_F*rWLS(nb).H(:,n);
        %Update velocity
        rWLS(nb).currVel = rWLS(nb).currVel + rWLS(nb).currAcc * dt;
        %Update position
        rWLS(nb).currPos = rWLS(nb).currPos + rWLS(nb).currVel * dt;
    end
end
%% plot
X = (1:length(rWLS(nb).Data(1,:)))*dt;
Y = zeros(length(X),1);
for n = 1:length(X)
    Y(n) = rWLS(nb).Data(1,n);
end

end

