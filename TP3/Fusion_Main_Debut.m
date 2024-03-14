close all;
%% Load data
Load_DATA 
%N = length(pos_GT);
%GPS_Data(:,2) = pos_GT(:,4) + normrnd(0,1,N,1);
%GPS_Data(:,3) = pos_GT(:,5) + normrnd(0,1,N,1);
%% plot in map
figure
% adjuste to figure size
h = 1000 ;w = 1900 ;set(gcf, 'Position', [0 0 w h]) ;movegui(gcf, 'center') 
Im_Saclay = imread('Saclay.png');
image(bounds_enu(:,1)',bounds_enu(:,2)'+280,Im_Saclay)
hold on;

%**************************************************************************
%%1- Tracer les donnees representant la position fournie par le recepteur 
%et la verite terrain.
%% plot data in (East,North)
plot(pos_GT(:,4),-pos_GT(:,5), 'Color', 'blue');
plot(GPS_Data(:,2),-GPS_Data(:,3), 'Color', 'cyan'); 
title('\fontsize{16}\color{blue}Ground Truth \color{cyan}GPS Data \color{Green}KF Estimation \color{red}KF Predictions');

%Generation donnees simulees a partir de la verite terrain
%**************************************************************************
%% 2-Initialiser positions, matrices de covariance ...
N_begin = 1;
N = length(pos_GT);
state_predict_vec = zeros(N,4);
a_rdn = 0.5;
% State : (X V_x Y V_y)
dt = GPS_Data(2,1) - GPS_Data(1,1);
State_old_p = [GPS_Data(2,2) ; ...   %Remplacer par X
              (GPS_Data(2,2) - GPS_Data(1,2))/dt;...    %Remplacer par V_x
              GPS_Data(2,3) ; ...   %Remplacer par Y
              (GPS_Data(2,3) - GPS_Data(1,3))/dt];       %Remplacer par V_y
N_state = length(State_old_p);
%2- Initialiser matrice de covariance du bruit de mesure
R = 4;   % A completer      % in meter                                                                 
Rk = (R^2)*eye(2) ;   % A completer      % in meter

% 2-Initialiser Matrice de covariance du vecteur d'etat, 
Pk_old = R^2 * eye(4,4) ;   % A completer      % in meter

%**************************************************************************
%7- Enregistrer les positions estimees et leur ecarts-types correspondants
KF_Estimation =[GPS_Data(N_begin,2)   GPS_Data(N_begin,3)    R  R;...      
                GPS_Data(N_begin+1,2) GPS_Data(N_begin+1,3)  R  R];         
last_k_valid=2;
state_predict_vec(1,:) = State_old_p';
state_predict_vec(1,:) = State_old_p';
%% Iteration principale
model=3;   % Modele d'innovation sur l'etat                                                                 

for k=N_begin+2:N
    V = sqrt(State_old_p(2)^2 + State_old_p(4)^2);
    theta = acos(State_old_p(2)/V);
    dt = GPS_Data(k,1) - GPS_Data(last_k_valid,1);
    last_k_valid=k;
    %3- Matrice de covariance du processus d'innovation. 
    Sigma_Acceleration = 5.0;     % in m/s   
    Sigma_theta = 10.0*pi/180;    % in rad/s    
    % Innovation process from acceleration [ax ay]--> State
    B  = [dt^2/2 0 ;dt 0 ;0 dt^2/2 ;0 dt] ;   
    % Innovation process from [v theta] --> acceleration [ax ay]
    B1 = [cos(theta) -V*sin(theta); sin(theta) V*cos(theta)];               
  switch model
      case 1       % Acceleration isotrope et constante         
    Q1 = [Sigma_Acceleration^2 0; 0 Sigma_Acceleration^2];
    Qk = B*Q1*B' ;
      case 2       % Acceleration anisotrope en theta et constante en V
    Q2 = [Sigma_Acceleration^2 0; 0 Sigma_theta^2];
    Q1 = B1*Q2*B1';
    Qk = B*Q1*B';
      case 3        % Acceleration anisotrope en theta et variable en 1/V
    Vseuil=3;     % En m/s
      if V>Vseuil % in rad/s en 1 sur V a partir de Vseuil
      Sigma_theta=Sigma_theta*Vseuil/V;  
      end
    Q2 = [Sigma_Acceleration^2 0; 0 Sigma_theta^2]; 
    Q1 = B1*Q2*B1'; 
    Qk = B*Q1*B';
  end
%3- Matrice de transition d'etat   
    Ak = eye(4);
    Ak(1,2) = dt;
    Ak(3,4) = dt;
    
%**3- Programmer l etape de prediction du filtre de Kalman.
  State_Predict =  Ak*State_old_p + B*normrnd(0,a_rdn,2,1);  % A completer
  Pk_predict = Ak*Pk_old*Ak' + Qk;   % A completer   
  
%**4- Definir la matrice d'observation. 
  Hk = [1 0 0 0;
        0 0 1 0];   % A completer      
  Zk =  Hk*State_Predict;   % A completer
  Sk = Hk*Pk_predict*Hk' + Rk; % A completer
   
%**5- Programmer l'etape d'estimation.  
  nuk = GPS_Data(k,2:3)' - Zk;
  Kk = Pk_predict * Hk' / Sk ;   % A completer    
  State_Estim = State_Predict + Kk*nuk ;   % A completer               
  Pk_new = (eye(4) - Kk*Hk) * Pk_predict * (eye(4) - Kk*Hk)';   % A completer      

%**6 Update Matrix and state for next iteration
  Pk_old = Pk_new ;   % A completer    
  State_old_p = State_Estim ;   % A completer    

%**7 Realiser divers Traces deroulants
  mark = 'cs'; %='ko'
  fen =1;
  thik = 1;                                 %Epaisseur du trace                                
  plot(pos_GT(k,4),-pos_GT(k,5),'bs','linewidth',thik); 
  plot(GPS_Data(k,2),-GPS_Data(k,3),mark,'linewidth',thik);
  plot(State_Predict(1),-State_Predict(3),'rs','linewidth',thik); 
  plot(State_Estim(1),-State_Estim(3),'gs');
%pause();
  Tracer=mod(k,5)==1;   %Tracer un ellipsoide tous les 5 iterations
  if Tracer   
  Fleche = [0 1 0.9 0.9 1; 0 0 -0.1 0.1 0]; %Trace d'une fleche                           
  beta = linspace(0,2*pi);                  %Trace d'une ellipsoide
  Cercle = [cos(beta); sin(beta)];
%Ellipsoide d'incertitude sur les data GPS
    plot(GPS_Data(k,2)+R*cos(beta),-GPS_Data(k,3)+R*sin(beta),'c');
%Tracer la vitesse predite
    V_P = [State_Predict(2) State_Predict(4);...
        -State_Predict(4) State_Predict(2)]*Fleche;
    plot(State_Predict(1)+V_P(1,:),-State_Predict(3)+V_P(2,:),'r-'); 
%Ellipsoide d'incertitude sur les positions predites
    [V,D] = eig(Sk);  %Valeurs et vecteurs propres
    Ellipse = fen*V*sqrt(D)*Cercle ;                                            
   plot(State_Predict(1)+Ellipse(1,:),-State_Predict(3)+Ellipse(2,:),'r-'); 
 %trace de la vitesse  estimee
    V_S = [State_Estim(2) State_Estim(4);...
        -State_Estim(4) State_Estim(2)]*Fleche*dt;
    plot(State_Estim(1)+V_S(1,:),-State_Estim(3)+V_S(2,:),'g-');
    [V,D] = eig(Hk * Pk_new * Hk'); %Extraire valeurs et vecteurs propres
    Ellipse = V*sqrt(D)*Cercle ;    %Incertitude sur positions predites
    plot(State_Estim(1)+Ellipse(1,:),-State_Estim(3)+Ellipse(2,:),'g-');
%Incertitude sur Vitesse  
    Hv = [0 1 0 0 ;0 0 0 1 ];   %Extraire vitesse
    [V,D] = eig(Hv * Pk_new  * Hv');  % -> l'ellipse Vitesse
    Ellipse = V*sqrt(D)*Cercle ;  
    plot(State_Estim(1)+State_Estim(2)*dt+Ellipse(1,:),...
    -State_Estim(3)-State_Estim(4)*dt+Ellipse(2,:),'g.');   
  end
  
%6- Trace d'un lien entre les variables de la meme iteration ...
% pour verification visuelle de l'association
  plot([pos_GT(k,4),GPS_Data(k,2),State_Estim(1),State_Predict(1)],...
    [-pos_GT(k,5),-GPS_Data(k,3),-State_Estim(3),-State_Predict(3)],'k-');
 if (pos_GT(k,7)-GPS_Data(k,1)~= 0)
     pos_GT(k,7)-GPS_Data(k,1)
 end
%*7- Enregistrer les estimees et leur ecarts-types correspondants, 
  KF_Estimation = [KF_Estimation ; ...
      State_Estim(1) State_Estim(3) sqrt(Pk_new(1,1)) sqrt(Pk_new(3,3))];
%  plot([KF_Estimation(i-1,1),KF_Estimation(i,1)],...
%      [-KF_Estimation(i-1,2),-KF_Estimation(i,2)],'g-');

  state_predict_vec(k,:) = State_Estim';
  pause(0.01);
end

plot(state_predict_vec(:,1),-state_predict_vec(:,3),'g'); 

%**************************************************************************
%7- Estimation statistique du data set 
%7- Verite terrain - estimee .
Ecart=[pos_GT(:,4)-KF_Estimation(:,1),-pos_GT(:,5)+KF_Estimation(:,2)];
Mean= mean(Ecart)
GroundTruthMinusEstim=std(Ecart)
%**7- Verite terrain - donnees GPS
Ecart=[pos_GT(:,4)-GPS_Data(:,2),-pos_GT(:,5)+GPS_Data(:,3)];               
moyenne  = mean(Ecart);        
Deviation_STD = std(Ecart)    
Deviation_STD2 = std2(Ecart)  
Norm= sqrt(norm(Ecart)/length(Ecart))
%**************************************************************************


