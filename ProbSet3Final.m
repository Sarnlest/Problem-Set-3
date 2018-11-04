%Peter Ferguson, EE209AS, 10/31/2018
clear all
close all
clc

syms Wr Wl Px Py Th k; %Right wheel angular velocity, left wheel angular velocity, robot x position, robot y position, theta
%%
%Generate actual and estimated starting positions
XAct=[rand()*0.2+0.15;rand()*0.3+0.225;rand()*2*pi]; %Initialize actual initial state somewhere near center
XActHistory=XAct;
XEst=XAct; %Initialize estimation with exact knowledge of actual position.
XEst=[rand()*0.5;rand()*0.75;rand*2*pi]; %Initialize randomly.  Comment out to calculate with known initial state
XEstHistory=XEst;
%%
%Initializing constants and matrices
WB=0.085; %Wheel base=0.085m
WRad=0.02; %Wheel radius=0.02m
DelT=0.033; %Delta T, sampling rate
A=eye(3); %A matrix.  The state does not change without input.
B=[DelT*WRad*cos(Th)/2, DelT*WRad*cos(Th)/2; DelT*WRad*sin(Th)/2, DelT*WRad*sin(Th)/2; DelT*WRad/WB, -DelT*WRad/WB]; %B Matrix.  Transformation matrix from U to X
f=FCalc(Px,Py,Th,Wr,Wl,A,B); %calculate the generic form of f
F=[diff(f(1),Px),diff(f(1),Py),diff(f(1),Th);diff(f(2),Px),diff(f(2),Py),diff(f(2),Th);diff(f(3),Px),diff(f(3),Py),diff(f(3),Th)]; %Calculate the generic form of F
P=1*eye(3); %initialize covariance at 1
%X=[Px;Py;Th]
%U=[Wr;Wl]
%Wall1=[k,0];% in meters for all k
%Wall2=[k,0.75];% in meters for all k
%Wall3=[0,k];% in meters for all k
%Wall4=[0.5,k];% in meters for all k

HEstHist=zeros(2,100);
%%
%Simulate
%Initialize U
IWr=randn()*50; %start with random velocity
IWl=randn()*50;
for l=1:100
    
    
    %%
    %Calculate actual state
    WSlip=randn(2,1)*(0.05*130*2*pi/60)^2; %wheel slip is 5% stdev normal noise on max wheel angular velocity (130rpm)
    XAct=subs(FCalc(XAct(1),XAct(2),XAct(3),IWr,IWl,A,B),Th,XAct(3))+subs(B,XAct(3))*WSlip; %Calculate the actual state
    XAct=double(XAct); %Change from symbolic to improve speed
    %%
    %Make sure it doesn't run out of bounds
    XAct(1)=min(XAct(1),0.499999);
    XAct(1)=max(XAct(1),0.000001);
    XAct(2)=min(XAct(2),0.7499999);
    XAct(2)=max(XAct(2),0.000001);
    
    %%
    %Calculate h for the actual state and corresponding actual measurements
    h1=Inf; %Initialize
    h2=Inf; %Initialize
    h1new=FLM(XAct(1),XAct(2),XAct(3)); %Calculate the forward distance to each wall
    h2new=RLM(XAct(1),XAct(2),XAct(3)); %Calculate the right distance to each wall
    for i=1:4 %Set h1 and h2 to be the forward and right measurements to the nearest wall in the correct direction.  Set H1Eq and H2Eq to the corresponding index number
        if double(h1new(i))>=0 && double(h1new(i))<h1
            h1=double(h1new(i));
            H1Eq=i;
        end
        if double(h2new(i))>=0 && double(h2new(i))<h2
            h2=double(h2new(i));
            H2Eq=i;
        end
    end
    h3=XAct(3);
    V=[h1*randn()*0.12;h2*randn()*0.12;randn()*3*pi/180];%Add sensor noise of stdev of 12% of current value for laser range finders, and 0.1degree/s for the angle based on the data sheets
    YAct=[h1;h2;h3]+V;
    
    %%Begin Kalman Filter for state estimation
    %%
    %Calculate time update state estimate and covariance based on process
    XEstNew=subs(FCalc(XEst(1),XEst(2),XEst(3),IWr,IWl,A,B),Th,XEst(3)); %Update state estimate based on previous state, inputs, and f
    F_t=subs(F,[Th,Wr,Wl],[XEstNew(3),IWr,IWl]); %Evaluate F at x (it only depends on Theta and U)
    Q_t=subs(B,XEstNew(3))*eye(2)*(0.05*130*2*pi/60)^2*subs(B,XEstNew(3))'; %Covariance of process noise
    PNext=F_t*P*F_t'+Q_t; %Calculate P_(t|t-1)
    
    %%
    %Calculate h and H
    h1=Inf; %Initialize
    h2=Inf; %Initialize
    h1new=FLM(XEstNew(1),XEstNew(2),XEstNew(3)); %Calculate the forward distance to each wall
    h2new=RLM(XEstNew(1),XEstNew(2),XEstNew(3)); %Calculate the right distance to each wall
    for i=1:4 %Set h1 and h2 to be the forward and right measurements to the nearest wall in the correct direction.  Set H1Eq and H2Eq to the corresponding index number
        if h1new(i)>=0 && h1new(i)<h1
            h1=h1new(i);
            H1Eq=i;
        end
        if h2new(i)>=0 && h2new(i)<h2
            h2=h2new(i);
            H2Eq=i;
        end
    end
    HEstHist(:,l)=[H1Eq;H2Eq];
    h3=XEstNew(3);
    YEst=[h1;h2;h3];
    H1=FLM(Px,Py,Th); %initialize the possible equations for H1
    H2=RLM(Px,Py,Th); %Initialize the possible equations for H2
    H1=H1(H1Eq); %Select the equation for the distance to the wall the forward laser actually hits
    H2=H2(H2Eq); %Select the equation for the distance to the wall the right laser actually hits
    H3=Th;
    H=[diff(H1,Px),diff(H1,Py),diff(H1,Th);diff(H2,Px),diff(H2,Py),diff(H2,Th);0,0,1]; %Find the partial derivatives of h wrt x
    H=H+[Px+Py,0,0;0,0,0;0,0,0]; %dummy addition to enable use of matlab function, removes later
    Ht=matlabFunction(H);
    H_t=Ht(XEstNew(1),XEstNew(2),XEstNew(3));
    H_t=H_t-[XEstNew(1)+XEstNew(2),0,0;0,0,0;0,0,0]; %eliminate dummy addition
    
    %%
    %Calculate difference between measured and expected Y, Kalman Gain
    YError=[YAct(1);YAct(2);YAct(3)]-YEst; %Difference between measured and expected Y
    R_t=[(YAct(1)*0.12)^2,0,0;0,(YAct(2)*0.12)^2,0;0,0,(3*pi/180)]; %assume a gray target outside around 40cm with 33ms pulses for 12% standard deviation on the laser ranges
    %No data on the datasheet for the IMU for variance or std dev of
    %magnetic north readings.  Used gyro noise of 0.1 degrees/second times DelT as standard dev.
    K_t=PNext*H_t'*(H_t*PNext*H_t'+R_t)^-1; %Kalman Gain
    %%
    %Calculate the estimated state with measurement considered.  Bound within
    %the rectangle
    XEst=vpa(XEstNew+K_t*YError,10);
    if XEst(1)>0.5
        XEst(1)=0.499999999999;
    elseif XEst(1)<0
        XEst(1)=0.000000000001;
    end
    if XEst(2)>0.75
        XEst(2)=0.7499999999999;
    elseif XEst(2)<0
        XEst(2)=0.000000000001;
    end
    
    %%
    %Update covariance
    P=vpa((eye(3)-K_t*H_t)*PNext*(eye(3)-K_t*H_t)'+K_t*R_t*K_t',10);
    %End Kalman Filter state estimate
    XError=double(XAct-XEst);
    DetP=det(P);
    XActHistory=cat(2,XActHistory,double(XAct));
    XEstHistory=cat(2,XEstHistory,double(XEst));
    
    
    %Update U randomly
    IWr=IWr+(5*randn());  %Randomly change input Wr 
    IWl=IWl+(5*randn());  %Randomly change input Wl 
    IWr=min(IWr,130*2*pi/60);  %MAX RPM of motor is 130.
    IWr=max(IWr,-130*2*pi/60);
    IWl=min(IWl,130*2*pi/60);
    IWl=max(IWl,-130*2*pi/60);
end
%%
%Plot Actual and Estimate
Time=size(XActHistory,2); %Get how many time steps
XActHistory(3,:)=mod(XActHistory(3,:),2*pi); %express theta between 0 and 2*pi
XEstHistory(3,:)=mod(XEstHistory(3,:),2*pi); %express theta between 0 and 2*pi
subplot(3,1,1)
plot(((1:Time)-1)*DelT,XActHistory(1,:)) %Plot actual X coordinate over time
hold on
plot(((1:Time)-1)*DelT,XEstHistory(1,:)) %Plot estimated X coordinate over time
ylabel('X (m)')
legend('Actual','Estimated')
axis([0,(Time-1)*DelT,0,0.5])
subplot(3,1,2)
plot(((1:Time)-1)*DelT,XActHistory(2,:)) %Plot actual Y coordinate over time
hold on
plot(((1:Time)-1)*DelT,XEstHistory(2,:)) %Plot estimated Y coordinate over time
legend('Actual','Estimated')
ylabel('Y (m)')
axis([0,(Time-1)*DelT,0,0.75])
subplot(3,1,3)
plot(((1:Time)-1)*DelT,XActHistory(3,:)) %Plot actual theta over time
hold on
plot(((1:Time)-1)*DelT,XEstHistory(3,:)) %Plot estimated theta over time
legend('Actual','Estimated')
ylabel('\theta (rad)')
xlabel('Time')
axis([0,(Time-1)*DelT,0,2*pi])
%%
%Plot which wall each laser is hitting.
figure
plot((2:Time)-1,HEstHist(1,:))
hold on
plot((2:Time)-1,HEstHist(2,:))
legend('H1','H2')

%%Functions
function [a]=FLM(Px,Py,Th)
a=[(0-Py)/sin(Th),(0.75-Py)/sin(Th),(0-Px)/cos(Th),(0.5-Px)/cos(Th)];%Distance of forward laser measurement to each wall 1-4 respectively
end
function [a]=RLM(Px,Py,Th)
a=[(0-Py)/sin(Th-pi/2),(0.75-Py)/sin(Th-pi/2),(0-Px)/cos(Th-pi/2),(0.5-Px)/cos(Th-pi/2)];%Distance of right laser measurement to each wall 1-4 respectively
end
function [a]=FCalc(Px,Py,Th,Wr,Wl,A,B)
a=A*[Px;Py;Th]+B*[Wr;Wl];
end