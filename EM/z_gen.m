 function zv = z_gen(N, M, Q)
 % zv = z_gen(N, M, Q)
 % Generates the Nth possible value of z(v), 0<=N<=M^Q-1
 %
 % INPUTS:
 % N - indexes range of z(v), 0<=N<=M^Q-1
 % M - number of MoG components
 % Q - number of ICs assumed to follow MoG distribution ('free' ICs)
 %
 % OUTPUT:
 % zv - Qx1 vector containing the (N+1)th possible value of z(v)

     if ( N < 0 || N >= M^Q )
         disp('error in input N: Must be between 0 and M^Q-1'); end;
     zv    = zeros(Q ,1);
     for q = 1:Q
           zv(q) =floor ( N / M^(Q-q) );
           N     = N - zv(q)*M^(Q-q);
     end;
     zv    = zv + 1;
 end