clear all;

DIFF_INC = 1e-4;

%pars = [D; B;                                   % Convex combination parameters 1-2
%    A;                                          % Wall strain parameter 3 
%    K_P; K_b; K_pb; K_pr; K_s;                  % Gains 4-8
%    tau_P; tau_b; taup_pb tau_pr; tau_s; tau_H; % Time Constants 9-14
%    q_w; q_pb; q_pr; q_s;                       % Sigmoid Steepnesses 15-18
%    xi_w; xi_pb; xi_pr; xi_s;                   % Sigmoid Shifts 19-22
%    HI; H_pb; H_pr; H_s;                        % Heart Rate Parameters 23-26
%    HIa; tchar];                                % HR mean after VM 27          

cd ../Sensitivities
N = 2;
for i = 1:N
   s = strcat('Ex',num2str(i),'_sens.mat');  
   load(s);
   
   INDMAP   = [11 19 23 24 25 27 28];
   sensN = sens(:,INDMAP);
   [U,S,V] = svd(sensN, 0);
   svals=diag(S);
   svals./svals(1);
   numopt = max(find(svals./svals(1) > DIFF_INC*10))
   [q,r,imap] = qr(V(:,1:numopt)',0);
   disp(strcat(num2str(i), 'subset:'));
   INDMAPn = sort(INDMAP(imap(1:numopt)))'
    
   S = sens(:,INDMAPn);
   [m,n] = size(S);

   A  = transpose(S)*S;
   Ai = inv(A);
   D  = diag(Ai);

   %disp('condition number of A = transpose(S)S and S');
   %disp([ cond(A) cond(S) ] );

   [a,b] = size(Ai);
   for i = 1:a
        for j = 1:b
            r(i,j)=Ai(i,j)/sqrt(Ai(i,i)*Ai(j,j)); % covariance matrix
        end;
   end;

   rn = triu(r,1); % extract upper triangular part of the matrix
   [i,j] = find(abs(rn)>0.90); % parameters with a value bigger than 0.85 are correlated
                                
   disp('correlated parameters');
   for k = 1:length(i)
       disp([INDMAPn(i(k)),INDMAPn(j(k)),rn(i(k),j(k))]);
   end
end;
 cd ../Core   