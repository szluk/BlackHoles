clear all

%% Black Hole Merger as an Event Converting Two Qubits Into One %%
%% BH Hamiltonians %%
% Based on
% https://www.researchgate.net/publication/391835509_Black_Hole_Merger_as_an_Event_Converting_Two_Qubits_Into_One
% (c) Szymon Lukaszyk
% licensed under MIT License
% email: szymon@patent.pl
% History
% v1: 17.05.2025
% v2: 26.07.2025 new version of eigenvectors

disp('%%% GENERAL QUBIT %%%')
syms o0 o1 o2 o3 real
H = [o0+o3   o1-i*o2;
     o1+i*o2 o0-o3];

% eigenvalues
o=sqrt(o1^2+o2^2+o3^2);
E0 = o0-o;
E1 = o0+o;

% eigenvectors
E0ket = [-(o1-i*o2); o3+o]/sqrt( 2*o*(o+o3) );
E1ket = [-(o1-i*o2); o3-o]/sqrt( 2*o*(o-o3) );
E0bra = [-(o1+i*o2) o3+o]/sqrt( 2*o*(o+o3) );
E1bra = [-(o1+i*o2) o3-o]/sqrt( 2*o*(o-o3) );
% alternative version congruent with BH Hamiltonian
E0ket = (o-o3)*[1; -(o1+i*o2)/(o-o3)]/sqrt( 2*o*(o-o3) );
E1ket = (o+o3)*[1;  (o1+i*o2)/(o+o3)]/sqrt( 2*o*(o+o3) );
E0bra = (o-o3)*[1  -(o1-i*o2)/(o-o3)]/sqrt( 2*o*(o-o3) );
E1bra = (o+o3)*[1   (o1-i*o2)/(o+o3)]/sqrt( 2*o*(o+o3) );

disp('% check 1 (orthonormality)')
chk=E0bra*E0ket;
chk=simplify(chk) % 1
chk=E0bra*E1ket;
chk=simplify(chk) % 0
chk=E1bra*E0ket;
chk=simplify(chk) % 0
chk=E1bra*E1ket;
chk=simplify(chk) % 1

disp('% check 2 (eigen equation)')
chk=H*E0ket - E0*E0ket;
chk=simplify(chk) % null vector
chk=H*E1ket - E1*E1ket;
chk=simplify(chk) % null vector

% check 3 (spectral decomposition)
chk=H - (E0*E0ket*E0bra + E1*E1ket*E1bra);
chk=simplify(chk) % null matrix

%return

disp('%%% 1st BH A %%%') 

syms A B th a b dta dtb dtab hbar real
% A=E_A, B=E_B       BH energies
% a=psi_A, b=psi_B   BH phases
% th=dt/hbar         simplifying factor

% orthogonalization interval
dto_A = hbar*pi/A;

% Hamiltonian
Ha = A*[1         exp(-i*a);
        exp(i*a)  1]/2;

%[Evec, Eval]=eig(Ha)

% Hamiltonian eigenvalues
E0 = 0;
E1 = A;

% Hamiltonian eigenvectors
E0ket = [1; -exp( i*a)]/sqrt(2);
E1ket = [1;  exp( i*a)]/sqrt(2);

E0bra = [1  -exp(-i*a)]/sqrt(2);
E1bra = [1   exp(-i*a)]/sqrt(2);

%other version
%E0ket = [-exp(-i*a); 1]/sqrt(2);
%E1ket = [ exp(-i*a); 1]/sqrt(2);

%E0bra = [-exp( i*a)  1]/sqrt(2);
%E1bra = [ exp( i*a)  1]/sqrt(2);

disp('% check 1 (orthonormality)')
E0bra*E0ket % 1
E0bra*E1ket % 0
E1bra*E0ket % 0
E1bra*E1ket % 1

disp('% check 2 (Schrodinger equation)')
Ha*E0ket - E0*E0ket % null vector
Ha*E1ket - E1*E1ket % null vector

disp('% check 3 (spectral decomposition)')
Ha - (E0*E0ket*E0bra + E1*E1ket*E1bra) % null matrix

disp('% unitary evolution operator of the Hamiltonian Ha')
%Ua = expm(-i*Ha*th)

%Ua1 = [ ( exp(-A*i*dta/hbar) + 1 )            ( (exp(-A*i*dta/hbar) - 1)*exp(-a*i) );
%       ( (exp(-A*i*dta/hbar) - 1)*exp(a*i) )  ( exp(-A*i*dta/hbar) + 1 )]/2;

Ua =  exp(-i*A*dta/(2*hbar)) * [  cos( A*dta/(2*hbar) )            -i*sin( A*dta/(2*hbar) )*exp(-a*i);
                                 -i*sin( A*dta/(2*hbar) )*exp(a*i)  cos( A*dta/(2*hbar) ) ];
%chk=Ua-Ua1
%chk=simplify(chk)

disp('% check 1 (unitarity)')
chk=Ua*Ua'
chk=simplify(chk) % identity matrix

disp('% check 2 (evolution at orthogonalization interval)') 
Uao = subs(Ua, dta, dto_A) % = subs(Ua, dta, -dto_A)
Uao2 = subs(Ua, dta, i*dto_A)
%Uao = [ 0         -exp(-a*i);
%       -exp(a*i)   0]

return

disp('%%% 2nd BH B %%%')
% orthogonalization interval
dto_B = hbar*pi/B;

% Hamiltonian
Hb = B*[1         exp(-i*b);
        exp(i*b)  1]/2;

% unitary evolution operator of the Hamiltonian Hb
%Ub = [ ( exp(-B*i*dtb/hbar) + 1 )             ( (exp(-B*i*dtb/hbar) - 1)*exp(-b*i) );
%       ( (exp(-B*i*dtb/hbar) - 1)*exp(b*i) )  ( exp(-B*i*dtb/hbar) + 1 )]/2;
Ub =  exp(-i*B*dtb/(2*hbar)) * [  cos( B*dtb/(2*hbar) )            -i*sin( B*dtb/(2*hbar) )*exp(-b*i);
                                 -i*sin( B*dtb/(2*hbar) )*exp(b*i)  cos( B*dtb/(2*hbar) ) ];
disp('% check commutation')
chk_a = kron(Ha, eye(2))
chk_b = kron(eye(2), Hb)
chk_a*chk_b - chk_b*chk_a % 0

disp('%%% A system of two independent non-interacting BHs A and B %%%')
% orthogonalization interval
dto_AB = hbar*pi/(A+B);

%Hab = kron(Ha, eye(2)) + kron(eye(2), Hb)

Hab = [A+B        B*exp(-i*b) A*exp(-i*a) 0;
       B*exp(i*b) A+B         0           A*exp(-i*a)
       A*exp(i*a) 0           A+B         B*exp(-i*b);
       0          A*exp(i*a)  B*exp(i*b)  A+B]/2 

%[Evec, Eval]=eig(Hab)

% eigenvalues
E0 = 0;
E1 = B;
E2 = A;
E3 = A+B;

% eigenvectors
E00ket = [1; -exp( i*b); -exp( i*a);  exp( i*(a+b))]/2;
E0Bket = [1;  exp( i*b); -exp( i*a); -exp( i*(a+b))]/2;
EA0ket = [1; -exp( i*b);  exp( i*a); -exp( i*(a+b))]/2;
EABket = [1;  exp( i*b);  exp( i*a);  exp( i*(a+b))]/2;

E00bra = [1  -exp(-i*b)  -exp(-i*a)   exp(-i*(a+b))]/2;
E0Bbra = [1   exp(-i*b)  -exp(-i*a)  -exp(-i*(a+b))]/2;
EA0bra = [1  -exp(-i*b)   exp(-i*a)  -exp(-i*(a+b))]/2;
EABbra = [1   exp(-i*b)   exp(-i*a)   exp(-i*(a+b))]/2;

disp('% check 1 (orthonormality)')
E00bra*E00ket % 1
E00bra*E0Bket % 0
EABbra*EABket % 1
EABbra*EA0ket % 0
%...

disp('% check 2 (Schrodinger equation)')
chk=Hab*E00ket - E0*E00ket; 
chk=simplify(chk) % null vector
chk=Hab*E0Bket - E1*E0Bket;
chk=simplify(chk) % null vector
chk=Hab*EA0ket - E2*EA0ket;
chk=simplify(chk) % null vector
chk=Hab*EABket - E3*EABket;
chk=simplify(chk) % null vector

disp('% check 3 (spectral decomposition)')
chk=Hab - (E0*E00ket*E00bra + E1*E0Bket*E0Bbra + E2*EA0ket*EA0bra + E3*EABket*EABbra)
chk=simplify(chk)  % null matrix

disp('% unitary evolution operator of the Hamiltonian Ha')
%Uab1 = expm(-i*Hab*th)
%Uab1=simplify(Uab1) 

% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

At2=(A*dtab/hbar)/2;
Bt2=(B*dtab/hbar)/2;

Uab =[( cos(At2)*cos(Bt2)*(cos(At2) - sin(At2)*i)*(cos(Bt2) - sin(Bt2)*i)                                         ) (-cos(At2)*sin(Bt2)*(cos(b) - sin(b)*i)*(cos(At2) - sin(At2)*i)*(cos(Bt2) - sin(Bt2)*i)*i                   ) (-cos(Bt2)*sin(At2)*(cos(a) - sin(a)*i)*(cos(At2) - sin(At2)*i)*(cos(Bt2) - sin(Bt2)*i)*i                   ) (-sin(At2)*sin(Bt2)*(cos(a) - sin(a)*i)*(cos(b) - sin(b)*i)*(cos(At2) - sin(At2)*i)*(cos(Bt2) - sin(Bt2)*i) );
      (-cos(At2)*sin(Bt2)*(cos(b) + sin(b)*i)*(cos(At2) - sin(At2)*i)*(cos(Bt2) - sin(Bt2)*i)*i                   ) ( cos(At2)*cos(Bt2)*(cos(At2) - sin(At2)*i)*(cos(Bt2) - sin(Bt2)*i)                                         ) (-sin(At2)*sin(Bt2)*(cos(a) - sin(a)*i)*(cos(b) + sin(b)*i)*(cos(At2) - sin(At2)*i)*(cos(Bt2) - sin(Bt2)*i) ) (-cos(Bt2)*sin(At2)*(cos(a) - sin(a)*i)*(cos(At2) - sin(At2)*i)*(cos(Bt2) - sin(Bt2)*i)*i                   );
      (-cos(Bt2)*sin(At2)*(cos(a) + sin(a)*i)*(cos(At2) - sin(At2)*i)*(cos(Bt2) - sin(Bt2)*i)*i                   ) (-sin(At2)*sin(Bt2)*(cos(a) + sin(a)*i)*(cos(b) - sin(b)*i)*(cos(At2) - sin(At2)*i)*(cos(Bt2) - sin(Bt2)*i) ) ( cos(At2)*cos(Bt2)*(cos(At2) - sin(At2)*i)*(cos(Bt2) - sin(Bt2)*i)                                         ) (-cos(At2)*sin(Bt2)*(cos(b) - sin(b)*i)*(cos(At2) - sin(At2)*i)*(cos(Bt2) - sin(Bt2)*i)*i                   );
      (-sin(At2)*sin(Bt2)*(cos(a) + sin(a)*i)*(cos(b) + sin(b)*i)*(cos(At2) - sin(At2)*i)*(cos(Bt2) - sin(Bt2)*i) ) (-cos(Bt2)*sin(At2)*(cos(a) + sin(a)*i)*(cos(At2) - sin(At2)*i)*(cos(Bt2) - sin(Bt2)*i)*i                   ) (-cos(At2)*sin(Bt2)*(cos(b) + sin(b)*i)*(cos(At2) - sin(At2)*i)*(cos(Bt2) - sin(Bt2)*i)*i                   ) ( cos(At2)*cos(Bt2)*(cos(At2) - sin(At2)*i)*(cos(Bt2) - sin(Bt2)*i)                                         )]

disp('% check 1 (unitarity)')
%chk=Uab*Uab';
%chk=simplify(chk) % identity matrix

% check 2 (determinant)
%chk=det(Uab);
%chk=simplify(chk) % OK

disp('% check 3 (tensor product of the individual evolution operators)')
%Uabk = kron(Ua,Ub);
%Uabk = subs(Uabk, dta, dtab);
%Uabk = subs(Uabk, dtb, dtab);
%Uabk = simplify(Uabk)
%chk=Uab - Uabk
%chk=simplify(chk) % null matrix

Uab1 = exp( -i*(A+B)*dtab/(2*hbar) )*[
 (   cos(At2)*cos(Bt2)               ) (-i*cos(At2)*sin(Bt2)*exp(-i*b)     ) (-i*sin(At2)*cos(Bt2)*exp(-i*a)     ) (  -sin(At2)*sin(Bt2)*exp(-i*(a+b)) );
 (-i*cos(At2)*sin(Bt2)*exp( i*b)     ) (   cos(At2)*cos(Bt2)               ) (  -sin(At2)*sin(Bt2)*exp(-i*(a-b)) ) (-i*sin(At2)*cos(Bt2)*exp(-i*a)     );
 (-i*sin(At2)*cos(Bt2)*exp( i*a)     ) (  -sin(At2)*sin(Bt2)*exp( i*(a-b)) ) (   cos(At2)*cos(Bt2)               ) (-i*cos(At2)*sin(Bt2)*exp(-i*b)     );
 (  -sin(At2)*sin(Bt2)*exp( i*(a+b)) ) (-i*sin(At2)*cos(Bt2)*exp( i*a)     ) (-i*cos(At2)*sin(Bt2)*exp( i*b)     ) (   cos(At2)*cos(Bt2)               )];

chk=Uab - Uab1
chk=simplify(chk) % null matrix
I=Uab1*Uab1';
I=simplify(I);

%{
At2=(A*dtab_ort/hbar)/2 = pi/2;
Bt2=(B*dtab_ort/hbar)/2 = pi/2;
cos(At2)=cos(Bt2)=0
sin(At2)=sin(Bt2)=1

Uab1 = [
0             0             0             exp(-i*(a+b));
0             0             exp(-i*(a-b)) 0;
0             exp( i*(a-b)) 0             0;
exp( i*(a+b)) 0             0             0];

Uab1_12 = [ Uab1(1,1)    Uab1(1,2);
            Uab1(2,1)    Uab1(2,2)];
I_12=Uab1_12*Uab1_12';
I_12=simplify(I_12)
I_12 =
[ cos((A*dtab)/(2*hbar))^2,                        0]
[                        0, cos((A*dtab)/(2*hbar))^2]
 
Uab1_34 = [ Uab1(3,3)    Uab1(3,4);
            Uab1(4,3)    Uab1(4,4)];
I_34=Uab1_34*Uab1_34';
I_34=simplify(I_34)
I_34 =
[ cos((A*dtab)/(2*hbar))^2,                        0]
[                        0, cos((A*dtab)/(2*hbar))^2]

Uab1_13 = [ Uab1(1,1)    Uab1(1,3);
            Uab1(3,1)    Uab1(3,3)];
I_13=Uab1_13*Uab1_13';
I_13=simplify(I_13)
I_13 =
[ cos((B*dtab)/(2*hbar))^2,                        0]
[                        0, cos((B*dtab)/(2*hbar))^2]

Uab1_24 = [ Uab1(2,2)    Uab1(2,4);
            Uab1(4,2)    Uab1(4,4)];
I_24=Uab1_24*Uab1_24';
I_24=simplify(I_24)
I_24 =
[ cos((B*dtab)/(2*hbar))^2,                        0]
[                        0, cos((B*dtab)/(2*hbar))^2]

Uab1_14 = [ Uab1(1,1)    Uab1(1,4);
            Uab1(4,1)    Uab1(4,4)];
I_14=Uab1_14*Uab1_14';
I_14=simplify(I_14)
I_14 =
[                                                                            (cos((A*dtab)/hbar)*cos((B*dtab)/hbar))/2 + 1/2, -2*cos((A*dtab)/(2*hbar))*cos((B*dtab)/(2*hbar))*sin((A*dtab)/(2*hbar))*sin((B*dtab)/(2*hbar))*(cos(a + b) - sin(a + b)*i)]
[ -2*cos((A*dtab)/(2*hbar))*cos((B*dtab)/(2*hbar))*sin((A*dtab)/(2*hbar))*sin((B*dtab)/(2*hbar))*(cos(a + b) + sin(a + b)*i),                                                                            (cos((A*dtab)/hbar)*cos((B*dtab)/hbar))/2 + 1/2]
=
[ 1/2, 0]
[ 0, 1/2]

Uab1_23 = [ Uab1(2,2)    Uab1(2,3);
            Uab1(3,2)    Uab1(3,3)];
I_23=Uab1_23*Uab1_23';
I_23=simplify(I_23)
I_23 =
[                                                                            (cos((A*dtab)/hbar)*cos((B*dtab)/hbar))/2 + 1/2, -2*cos((A*dtab)/(2*hbar))*cos((B*dtab)/(2*hbar))*sin((A*dtab)/(2*hbar))*sin((B*dtab)/(2*hbar))*(cos(a - b) - sin(a - b)*i)]
[ -2*cos((A*dtab)/(2*hbar))*cos((B*dtab)/(2*hbar))*sin((A*dtab)/(2*hbar))*sin((B*dtab)/(2*hbar))*(cos(a - b) + sin(a - b)*i),                                                                            (cos((A*dtab)/hbar)*cos((B*dtab)/hbar))/2 + 1/2]

%}


return

%chk=Uab(1,1)-Uab1(1,1) %0
%chk=Uab(1,2)-Uab1(1,2) %0
%chk=Uab(1,3)-Uab1(1,3) %0
%chk=Uab(1,4)-Uab1(1,4) %0

%chk=Uab(2,1)-Uab1(2,1) %0
%chk=Uab(2,2)-Uab1(2,2) %0
%chk=Uab(2,3)-Uab1(2,3) %0
%chk=Uab(2,4)-Uab1(2,4) %0

%chk=Uab(3,1)-Uab1(3,1) %0
%chk=Uab(3,2)-Uab1(3,2) %0
%chk=Uab(3,3)-Uab1(3,3) %0
%chk=Uab(3,4)-Uab1(3,4) %0

%chk=Uab(4,1)-Uab1(4,1) %0
%chk=Uab(4,2)-Uab1(4,2) %0
%chk=Uab(4,3)-Uab1(4,3) %0
%chk=Uab(4,4)-Uab1(4,4) %0
%chk=simplify(chk)

disp('% Uab at Orthogonalization intervals')
Uabo = [ 0             0             0             exp(-i*(a+b));
         0             0             exp(-i*(a-b)) 0; 
         0             exp( i*(a-b)) 0             0;
         exp( i*(a+b)) 0             0             0];

disp('% check 1 (unitarity)')
Uabo*Uabo'
det(Uabo)

return
