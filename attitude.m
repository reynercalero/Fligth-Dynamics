function [t,phi,theta,psi,phi_dt,theta_dt,psi_dt]=attitude();

data = dlmread('tello_imu_example.csv',',',1,0);

t = data(:,1);
p = data(:,2);
q = data(:,3);
r = data(:,4);
N = length(t);

phi = zeros(N,1);
theta = zeros(N,1);
psi = zeros(N,1);

phi_dt=zeros(N,1);
theta_dt=zeros(N,1);
psi_dt=zeros(N,1);


for i = 2:N;
   dt = t(i) - t(i-1);
   ph = phi(i-1);
   th = theta(i-1);
   if abs(cos(th)) < 1e-6;
     th=th + 0.001;

   endif

   H = [1,  sin(ph)*tan(th),  cos(ph)*tan(th);
         0,  cos(ph),         -sin(ph);
         0,  sin(ph)/cos(th),  cos(ph)/cos(th)];

   body_rates = [p(i),q(i),r(i)]';

   euler_rates = H * body_rates;

   phi_dt(i) = euler_rates(1);
   theta_dt(i) = euler_rates(2);
   psi_dt(i) = euler_rates(3);

   phi(i) = phi(i-1) + phi_dt(i)*dt;
   theta(i) = theta(i-1) + theta_dt(i)*dt;
   psi(i) = psi(i-1) + psi_dt(i)*dt;
 end


