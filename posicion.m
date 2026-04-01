function [t,px,py,pz,Vx,Vy,Vz,Xac,Yac,Zac,a_ned,a_real]=posicion();
  [t,phi,theta,psi]=attitude();
  data = dlmread('tello_imu_example.csv',',',1,0);

  Xac = data(:,5);
  Yac = data(:,6);
  Zac = data(:,7);
  N = length(t);

  Vx=zeros(N,1);
  Vy=zeros(N,1);
  Vz=zeros(N,1);

  px=zeros(N,1);
  py=zeros(N,1);
  pz=zeros(N,1);

  g=[0;0;-9.81];

  for i = 2:N;
    dt=t(i)-t(i-1);
    R=matriz_rotacion(phi(i),theta(i),psi(i));

    a_b=[Xac(i);Yac(i);Zac(i)];
    a_ned = R'*a_b;
    a_real=a_ned-g;

    Vx(i) = Vx(i-1) + a_real(1) * dt;
    Vy(i) = Vy(i-1) + a_real(2) * dt;
    Vz(i) = -(Vz(i-1) + a_real(3) * dt);

    px(i) = px(i-1) + Vx(i) * dt;
    py(i) = py(i-1) + Vy(i) * dt;
    pz(i) = pz(i-1) + Vz(i) * dt;
  end
end

