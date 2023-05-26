%CODIGO PARA IMPLEMENTAR EL CONTROLADOR
%ALUMNO: HERMAN RAYMUNDO SALVATIERRA 
%TP 2 item 1
%Profesor Pucheta
clear all, close all, clc
%Variables
Laa= 5e-3;
J= 0.004;
Ra= 0.2;
B= 0.005; %0,005
Ki= 6.5e-5;
Km= 0.055;
%Torque (Dato)
T1=1.15e-3; %Va a ir variando para pi/2 y -pi/2

%Planteo matrices A B C D de estado a partir de las ecuaciones del modelo
%siendo
% x1 = ia    :la corriente que circula
% x2 = wr    :la velocidad angular
% x3 = theta :el angulo posicion del motor
Mat_A = [-Ra/Laa -Km/Laa 0; Ki/J -B/J 0; 0 1 0];
Mat_B = [1/Laa; 0; 0];
Mat_C = [0 0 1];
Mat_D = [0];
%Matrices ampliadas
Mat_Aa = [Mat_A zeros(3,1); -Mat_C 0];
Mat_Ba = [Mat_B; 0];
Mat_Cc = [Mat_C 0];
%Se hacen varios intentos para el controlador por asignacion de polos pero
%me traen muchos problemas, y por eso decidí implementar el controlador por
%LQR por ser mas eficiente
%Despues de probar con distintos valores opté por el que me funciona bien

Q = diag([1 1/10000 100 10000000000000]);
 
R=100.0; 

%Construccion del Hamiltoniano para el calculo del controlador
Ha=[Mat_Aa -Mat_Ba*inv(R)*Mat_Ba'; -Q -Mat_Aa'];
[n,va]=size(Ha); % estas se definen y las usaré como limites para el for

[V,D]=eig(Ha);% son los 2 polos de la funcion de transf 
              %a lazo abierto del hamiltoniano
MX1X2=[]; %empieza vacia y luego recorre y se van guardando
          %los valores si es que tienen parte real negativa
for ii=1:n
    if real(D(ii,ii))<0 %recorro la diagonal principal buscando polos estables
        MX1X2=[MX1X2 V(:,ii)];
    end
end

MX1=MX1X2(1:n/2,:); % A esta variable MX1 le asigno el valor de la posicion n/2
MX2=MX1X2(n/2+1:end,:);% A esta variable MX2 le asigno el valor de la posicion n/2+1

P=real(MX2*inv(MX1));%Tomo la parte real por un tema de residuos
Ka=inv(R)*Mat_Ba'*P;
K=Ka(1:3); 
KI=-Ka(4);
eig(Mat_Aa-Mat_Ba*Ka)%Verifico polos con parte real negativa
% break
%Fin cálculo del controlador
J_(1)=0
V_(1)=0;
psi(1)=0;

%Para la simulación
delta_t=1e-3; %1e-5     ESTE PARAMETRO HACE QUE LA CORRIENTE CAMBIE

tiempo=10;
pasos=round(tiempo/delta_t);
Ci=[0 0 0 0];
t=0:delta_t:(tiempo-delta_t);

%Referencia
%cambie el 0.6 por 4
ref=(pi/2)*square(2*pi*t/4);

%Torque (original 0.6)
%probe cambiar esta formula (debo aumentar el 0.6)
T11=(T1/2)+(T1/2)*square(2*pi*t/4 );%Función para ir variando el torque
%T1=(1.15e-3/2)+(1.15e-3/2)*square(2*pi*t/0.6); 2.4335e-3

%Condiciones iniciales
x=zeros(4,pasos);
x(1,1)=Ci(1);%ia
x(2,1)=Ci(2);%w
x(3,1)=Ci(3);%theta
x(4,1)=Ci(4);
ua(1)=0;

for i=2:1:pasos
    estado=[x(1,i-1);x(2,i-1);x(3,i-1);x(4,i-1)];%Guardo el estado
    integracion=x(4,i-1)+delta_t*(ref(1,i-1)-Mat_Cc*estado);%Es el termino de la pseudo inversa
    u_actual=-Ka*estado;%(1:3)+integracion*KI;
    ua=[ua u_actual];
    
    %Ecuaciones del motor
    ia_p=(-Ra/Laa)*estado(1)-(Km/Laa)*estado(2)+(1/Laa)*u_actual;
    w_p=(Ki/J)*estado(1)-(B/J)*estado(2);%-(1/J)*T11(i-1);
    theta_p=estado(2);
    
    xp_actual=[ia_p; w_p; theta_p];
    
    xsig=estado(1:3)+delta_t*xp_actual;
    x(1,i)=xsig(1);
    x(2,i)=xsig(2);
    x(3,i)=xsig(3);
    x(4,i)=integracion;
end

figure
plot(t,x(3,1:end));title('Angulo Theta y referencia');hold on;grid on;
plot(t,ref);legend('Angulo Theta','ref');
%break
figure
plot(t,x(2,:));title('Velocidad angular, w');hold on;

figure
plot(t,x(1,:));title('Corriente Ia');hold on;

figure
plot(t,ua);title('Accion de control');hold on;

figure
plot(t,T11);title('Torque');hold on;
