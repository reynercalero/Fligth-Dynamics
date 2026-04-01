function [alpha, beta, gamma,t] =angulos_aeros()
    [t,px,py,pz,Vx,Vy,Vz,Xac,Yac,Zac,a_ned,a_real]=posicion();
    [t,phi,theta,psi,phi_dt,theta_dt,psi_dt]=attitude();


N = length(t);
alpha=zeros(N,1);
beta=zeros(N,1);
gamma=zeros(N,1);

for i=1:N
  R = matriz_rotacion(phi(i), theta(i), psi(i));
  Vned=[Vx(i);Vy(i);Vz(i)];
  Vb = R*Vned;

  u = Vb(1); v=Vb(2); w=Vb(3);
  Vtotal=norm(Vb);

  alpha(i)=rad2deg(atan2(-w,u));

  if Vtotal < 1e-6
     beta(i)=0;

  else
     beta(i)=rad2deg(asin(v/Vtotal));
  end

  gamma(i)=rad2deg(atan2(w,(u^2+v^2)^0.5));
endfor

end

