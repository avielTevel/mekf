%temp_velKF
% constants
close all

dT = 1;
T=11;
% set up desired trajectory

t = [0:dT:9]';
LL=length(dT);
posd = 1-cos(pi*t);
veld = t*0+1;
accd = pi^2*cos(pi*t);
postrue = posd*0;
veltrue = posd*0;
measurement = posd*0;
torque = posd*0;
pose = posd*0;
vele = posd*0;
posd=[-2 -1 0 1 2 3 4 5 6 7]'

% for ii = 1:1:LL
% posd(ii) = 0;
% veld(ii) = 0;
% accd(ii) = 0;
% end


% Show how differentiation blows up noise

mpos = posd + 0.01*randn(size(posd));
plot(t,mpos)
title('position with added noise')
xlabel('sample (100Hz)')


% plot(t,diffxy(t,mpos))                  
hold
plot(t,veld,'r')
hold
legend('differentiated position','true velocity')
xlabel('sample (100Hz)')
% Create idealized system to match filter assumptions
% state is (position, velocity)'
% We are assuming Q, R, and P0 to be diagonal

A = [ 1 1
      0 1 ]

B = [ 0
      0 ]

C = [ 1 0 ]

% process variance
Q = [ 1 0
      0 1]
% sensor nose variance
R = [ 1]

% initial state estimate variance
P0 = [ 5 0
       0 5]

% % % % Create some data
% % % state = [ sqrt(P0(1,1))*randn(1)
% % %           sqrt(P0(2,2))*randn(1) ];
% % % for ii = 1:length(posd)
% % %  postrue(ii) = state(1);
% % %  veltrue(ii) = state(2);
% % %  % simulate noisy measurement
% % %  measurement(ii) = C*state + sqrt(R(1,1))*randn(1);
% % %  torque(ii) = accd(ii);
% % %  process_noise = [ sqrt(Q(1,1))*randn(1)
% % %                    sqrt(Q(2,2))*randn(1) ];
% % %  state = A*state + B*accd(ii) + process_noise;
% % % end

postrue = posd;
veltrue = veld;
measurement = mpos;
torque = accd;



% Design filter
% Note that we can design filter in advance of seeing the data.
% % % Pm = P0;
% % % for ii = 1:1000
% % %  % measurement step
% % %  S = C*Pm*C' + R;
% % %  K = Pm*C'*inv(S);
% % %  Pp = Pm - K*C*Pm;
% % %  % prediction step
% % %  Pm = A*Pp*A' + Q;
% % % end
% % % 
% % % % Run the filter to create example output
% % % sem = [ 0 0 ]'
% % % for ii = 1:length(posd)
% % %  % measurement step
% % %  sep = sem + K*(measurement(ii)-C*sem);
% % %  pose(ii) = sep(1);
% % %  vele(ii) = sep(2);	   
% % %  % prediction step
% % %  sem = A*sep + B*torque(ii);
% % % end

a=zeros(2,10)
b=zeros(1,10)
c=zeros(1,10)
d=zeros(1,10)
e=zeros(1,10)

Pm = P0;
for ii = 1:10
 % measurement step
 S = C*Pm*C' + R;
 K = Pm*C'*inv(S);
 Pp = Pm - K*C*Pm;
 % prediction step
 Pm = A*Pp*A' + Q;
 a(:,ii) = K;
 b(ii)=Pm(1,1);
 c(ii)=Pm(1,2);
 d(ii)=Pm(2,1);
 e(ii)=Pm(2,2);

 disp(num2str(Pm'))
 disp('')
end

% Run the filter to create example output
sem = [ 0 0 ]'
for ii = 1:length(posd)
 % measurement step
 sep = sem + K*(measurement(ii)-C*sem);
 pose(ii) = sep(1);
 vele(ii) = sep(2);	   
 % prediction step
 sem = A*sep + B*torque(ii);
end



% Let's plot the Kalman filter output
ii = 1:length(pose);
plot(t,pose,'b',t,postrue,'r')
hold on
% plot(t,measurement,'b',t,postrue,'r')
plot(t,vele,'g',t,veltrue,'k')
legend('KF position','true position','KF velocity','true velocity')
xlabel('sample (k)')
return
% Let's compare to directly filtering the output.
vel1 = diff( measurement )/dT;
% get length right
vel1(length(veltrue)) = 0;
plot(t,ii,vel1,'b-.',ii,vele,'r')

% How well can we do with a first order filter?
[B,A]=butter(1,0.1);
vel2 = filter(B,A,vel1);
%plot(t,ii,vel2,'b',ii,veltrue,'m',ii,vele,'r')
plot(t,ii,vel2,'b',ii,veltrue,'r')
legend('1st order filtered velocity','true velocity')
xlabel('sample (100Hz)')
