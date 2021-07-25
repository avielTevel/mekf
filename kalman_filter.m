%temp_velKF
%constants
close all
%set up the time step to 0.1 seconds and the full simulation time to 6000
dT = 0.1; 
T=60;
%Numerical data
One_Three_Matrix=transpose([ 1 1 1]);
omega_t=One_Three_Matrix*sin(2*pi*t/150); %deg/sec
Epsilon_Eta_Square=10^-13;%rad^2/sec
Epsilon_n_Square=10^-11;%rad^2/sec^3
Epsilon_b=100; %arcsec
e1=0.378;
e2=-0.378;
e3= 0.756;
q=0.3780;
q0=transpose([e1 , e2 , e3 , q]);
Mu=3.1^-5*One_Three_Matrix; %rad/sec
r1=transpose([1 0 0]); %star tracker 1
r2=transpose([0 1 0]); %star tracker 2
%End of numerical data
%Task 1 - Run AEKF 1,2,3 and MEKF
X_gag_Initial=[0, 0, 0, 1, 0, 0, 0]; %X_gag_0/0
P_Initial=5*ones(7);
P_k_k=P_Initial;
Measurment=1;
for t=0:0.1:T
    %AEKF1
    %1)initialization the variables: q_gag_0/0 and P_0/0
    q_gag_k_k=X_gag_Initial;
    %Time promotion of the variables q_k/k and p_k/k
    %Defining the white noise matrices (The noise of the Gyro and the noise of the star trackers)
    Q_d_k=Epsilon_Eta_Square*dT*ones(3)+Epsilon_n_Square*dT*ones(3);
    ex=[0 -e3 e2
        e3 0  -e1
        -e2 e1 0];
    e=[e1 , e2, e3];
    Xi_qt=[ex+q*ones(3)
             -(e)];
    G_gag_k=-0.5*Xi_qt;
    %Defining phi where Omega is the triad gyro measurements
    Omega_t=omega_t(:,Measurment)+transpose(ones(1,3)*normrnd(0,Epsilon_n_Square)); %white noise with 0 mean and Epsilon n square intensity 
    phi=power(e,-0.5*Omega_t*dT); %dimentions?
    q_gag_kplus1_k=phi.*q_gag_k_k; 
    P_kplus1_kplus1=phi*P_k_k*transpose(phi)+G_gag_k*Q_d_k*tanspose(G_gag_k);
    %updating q_gag_k+1/k and P_k+1/k - updating by the star tracker
    b_kplus1=0;
    R_kplus1=0;
    H_gag_kplus1=H(q_gag_kplus1_k);
    K_kplus1=P_gag_kplus1_k*transpose(H_gag_kplus1)*(H_gag_kplus1*P_gag_kplus1_k*transpose(H_gag_kplus1)+R_kplus1)^-1;
    q_gag_kplus1_kplus1=q_gag_kplus1_k+K_kplus1(b_kplus1-A(q_gag_kplus1_k)*r_kplus1); %r is for a specific star , q_gag_kplus1_k size must be 1 to be ortogonal
    P_kplus1_kplus1=(I-K_kplus1*H_gag_kplus1)*P_gag_kplus1_k*transpose(I-K_kplus1*H_gag_kplus1)+K_kplus1*R_kplus1*transpose(K_kplus1);





    Measurment=Measurment+1;
end
% Show how differentiation blows up noise
% 
% mpos = posd + 0.01*randn(size(posd));
% plot(t,mpos)
% title('position with added noise')
% xlabel('sample (100Hz)')
% 
% 
% % plot(t,diffxy(t,mpos))                  
% hold
% plot(t,veld,'r')
% hold
% legend('differentiated position','true velocity')
% xlabel('sample (100Hz)')
% % Create idealized system to match filter assumptions
% % state is (position, velocity)'
% % We are assuming Q, R, and P0 to be diagonal
% 
% A = [ 1 1
%       0 1 ]
% 
% B = [ 0
%       0 ]
% 
% C = [ 1 0 ]
% 
% % process variance
% Q = [ 1 0
%       0 1]
% % sensor nose variance
% R = [ 1]
% 
% % initial state estimate variance
% P0 = [ 5 0
%        0 5]
% 
% % % % % Create some data
% % % % state = [ sqrt(P0(1,1))*randn(1)
% % % %           sqrt(P0(2,2))*randn(1) ];
% % % % for ii = 1:length(posd)
% % % %  postrue(ii) = state(1);
% % % %  veltrue(ii) = state(2);
% % % %  % simulate noisy measurement
% % % %  measurement(ii) = C*state + sqrt(R(1,1))*randn(1);
% % % %  torque(ii) = accd(ii);
% % % %  process_noise = [ sqrt(Q(1,1))*randn(1)
% % % %                    sqrt(Q(2,2))*randn(1) ];
% % % %  state = A*state + B*accd(ii) + process_noise;
% % % % end
% 
% postrue = posd;
% veltrue = veld;
% measurement = mpos;
% torque = accd;
% 
% 
% 
% % Design filter
% % Note that we can design filter in advance of seeing the data.
% % % % Pm = P0;
% % % % for ii = 1:1000
% % % %  % measurement step
% % % %  S = C*Pm*C' + R;
% % % %  K = Pm*C'*inv(S);
% % % %  Pp = Pm - K*C*Pm;
% % % %  % prediction step
% % % %  Pm = A*Pp*A' + Q;
% % % % end
% % % % 
% % % % % Run the filter to create example output
% % % % sem = [ 0 0 ]'
% % % % for ii = 1:length(posd)
% % % %  % measurement step
% % % %  sep = sem + K*(measurement(ii)-C*sem);
% % % %  pose(ii) = sep(1);
% % % %  vele(ii) = sep(2);	   
% % % %  % prediction step
% % % %  sem = A*sep + B*torque(ii);
% % % % end
% 
% a=zeros(2,10)
% b=zeros(1,10)
% c=zeros(1,10)
% d=zeros(1,10)
% e=zeros(1,10)
% 
% Pm = P0;
% for ii = 1:10
%  % measurement step
%  S = C*Pm*C' + R;
%  K = Pm*C'*inv(S);
%  Pp = Pm - K*C*Pm;
%  % prediction step
%  Pm = A*Pp*A' + Q;
%  a(:,ii) = K;
%  b(ii)=Pm(1,1);
%  c(ii)=Pm(1,2);
%  d(ii)=Pm(2,1);
%  e(ii)=Pm(2,2);
% 
%  disp(num2str(Pm'))
%  disp('')
% end
% 
% % Run the filter to create example output
% sem = [ 0 0 ]'
% for ii = 1:length(posd)
%  % measurement step
%  sep = sem + K*(measurement(ii)-C*sem);
%  pose(ii) = sep(1);
%  vele(ii) = sep(2);	   
%  % prediction step
%  sem = A*sep + B*torque(ii);
% end
% 
% 
% 
% % Let's plot the Kalman filter output
% ii = 1:length(pose);
% plot(t,pose,'b',t,postrue,'r')
% hold on
% % plot(t,measurement,'b',t,postrue,'r')
% plot(t,vele,'g',t,veltrue,'k')
% legend('KF position','true position','KF velocity','true velocity')
% xlabel('sample (k)')
% return
% % Let's compare to directly filtering the output.
% vel1 = diff( measurement )/dT;
% % get length right
% vel1(length(veltrue)) = 0;
% plot(t,ii,vel1,'b-.',ii,vele,'r')
% 
% % How well can we do with a first order filter?
% [B,A]=butter(1,0.1);
% vel2 = filter(B,A,vel1);
% %plot(t,ii,vel2,'b',ii,veltrue,'m',ii,vele,'r')
% plot(t,ii,vel2,'b',ii,veltrue,'r')
% legend('1st order filtered velocity','true velocity')
% xlabel('sample (100Hz)')
