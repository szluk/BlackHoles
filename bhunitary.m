clear all

%% Black Hole Merger as an Event Converting Two Qubits Into One %%
%% BH Unitary transformations %%
% Based on
% https://www.researchgate.net/publication/391835509_Black_Hole_Merger_as_an_Event_Converting_Two_Qubits_Into_One
% (c) Szymon Lukaszyk
% licensed under MIT License
% email: szymon@patent.pl
% History
% v1: 17.08.2025

% a    - phase factor
% A    - BH energy
% dta  - background time
% hbar - reduced Planck constant

syms a A dta hbar

% Hamiltonian eigenvectors
E0ket_b = [1; -exp( i*a)]/sqrt(2);
E1ket_b = [1;  exp( i*a)]/sqrt(2);

E0bra_b = [1  -exp(-i*a)]/sqrt(2);
E1bra_b = [1   exp(-i*a)]/sqrt(2);

%phi = A*dta/hbar;
syms phi
phio= pi;

% general unitary transformation
Ua =  exp(-i*phi/2) * [   cos( phi/2 )           -i*sin( phi/2 )*exp(-a*i);
                       -i*sin( phi/2 )*exp(a*i)     cos( phi/2 ) ];

% orthogonal unitary transformation with dta = dto_A = (hbar*pi/A)
Uao = -[ 0        exp(-a*i);
         exp(a*i) 0];

% Hamiltonian eigenvectors

ket0_b  = [1;  0];         % |0> in computational basis
ket1_b  = [0;  1];         % |1> in computational basis
ketneg_b= [1; -1]/sqrt(2); % |-> in computational basis
ketpos_b= [1;  1]/sqrt(2); % |+> in computational basis

% |E_0> eigenvector

disp('% check |E_0> eigenvector: should be null vectors')
E0ket= E0ket_b;
for n=1:5
    E0ket = Ua*E0ket;
    E0ket_n = E0ket; % for all dta, all n
    dif=simplify(E0ket-E0ket_n)
end

% |E_1> eigenvector

disp('% check |E_1> eigenvector: should be null vectors')
E1ket = E1ket_b;
for n=1:5
    E1ket = Ua*E1ket;
    E1ket_n = exp(-i*n*phi)*E1ket_b;
    dif=simplify(E1ket-E1ket_n)
end

disp('% check |E_1> eigenvector ortho interval: should be null vectors')
E1ket = E1ket_b;
for n=1:5
    E1ket = Uao*E1ket;
    E1ket_n = (-1)^n * E1ket_b;
    dif=simplify(E1ket-E1ket_n)
end

% |0> in computational basis

disp('% check |0> in computational basis: should be null vectors')
ket0 = ket0_b;
for n=1:5
    ket0 = Ua*ket0;
    ket0_n = [ exp(-i*n*phi)+1; exp(i*a) * (exp(-i*n*phi)-1) ]/2;
    dif=simplify(ket0-ket0_n)
end

disp('% check |0> in computational basis ortho interval: should be null vectors')
ket0 = ket0_b;
for n=1:5
    ket0 = Uao*ket0;
    if mod(n, 2) == 0 % n is even
        ket0_n = ket0_b;
    else % n is odd
        ket0_n = -exp(i*a)*ket1_b;
    end
    dif=simplify(ket0-ket0_n)
end

% |1> in computational basis

disp('% check |1> in computational basis: should be null vectors')
ket1 = ket1_b;
for n=1:5
    ket1 = Ua*ket1;
    ket1_n = [ exp(-i*a) * ( exp(-i*n*phi)-1); exp(-i*n*phi) + 1 ]/2;
    dif=simplify(ket1-ket1_n)
end

disp('% check |1> in computational basis ortho interval: should be null vectors')
ket1 = ket1_b;
for n=1:5
    ket1 = Uao*ket1;
    if mod(n, 2) == 0 % n is even
        ket1_n = ket1_b;    
    else % n is odd
        ket1_n = -exp(i*a)*ket0_b;
    end
    dif=simplify(ket0-ket0_n)
end

% |-> in computational basis

disp('% check |-> in computational basis: should be null vectors')
ketneg = ketneg_b;
for n=1:5
    ketneg = Ua*ketneg;
    ketneg_n = exp(-i*n*phi/2)*[ cos(n*phi/2) + i*sin(n*phi/2)*exp(-i*a); -cos(n*phi/2) - i*sin(n*phi/2)*exp(i*a)]/sqrt(2);    
    dif=simplify(ketneg-ketneg_n)
end

% |+> in computational basis

disp('% check |+> in computational basis: should be null vectors')
ketpos = ketpos_b;
for n=1:5
    ketpos = Ua*ketpos;
    ketpos_n =  exp(-i*n*phi/2)*[  cos(n*phi/2) - i*sin(n*phi/2)*exp(-i*a); cos(n*phi/2) - i*sin(n*phi/2)*exp(i*a)]/sqrt(2);
    dif=simplify(ketpos-ketpos_n)
end
