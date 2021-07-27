clear all

% Properties of pi-bit and 4*pi-bit black holes
%
% (c) Szymon £ukaszyk
% email: szymon@patent.pl
% licensed under MIT License.
% History
% 15.07.2021 1st version
% 25.07.2021 2nd version (corrected error in angle C calculation)
% 27.07.2021 cotan Laplacians

% pi-bit black hole
% Spherical triangles angles
A = 2*pi/3;
B = 2 + pi/6; 
C = 5*pi/3 - 4;
% C = 10*B-24

% Cosine rule for angles
cosa = (cos(B)^2-0.5)/sin(B)^2;
cosb = cot(B)/sqrt(3); sinb = sqrt(1-cosb^2);
cosc = (cos(C) + cos(C)^2)/sin(C)^2;
% cosa = cosc
% sin(B)^2 = sin(5*B - 12)^2

% Spherical edges (angles in radians)
a = acos(cosa); % a = 1.059395540392039, adeg = 60.698893299444990
b = acos(cosb); % b = 2.518636598528195, bdeg = 144.3072472228511
c = acos(cosc); % c = 1.059395540392039, cdeg = 60.698893299444990
% a = c

% Cartesian coordinates on unit sphere (r=1)
p1 = [ 0       0              1   ];
p2 = [ sinb    0              cosb];
p3 = [-sinb/2  sqrt(3)*sinb/2 cosb];
p4 = [-sinb/2 -sqrt(3)*sinb/2 cosb];
v=[p1;p2;p3;p4];

% Voronoi vertices
c123 = cross(p2-p1, p3-p1)/norm(cross(p2-p1, p3-p1));
c134 = cross(p3-p1, p4-p1)/norm(cross(p3-p1, p4-p1));
c214 = cross(p1-p2, p4-p2)/norm(cross(p1-p2, p4-p2));
c243 = cross(p4-p2, p3-p2)/norm(cross(p3-p2, p4-p2));
%c123 = ( 0.493644698009721,  0.855017697839831, 0.158932842759347)
%c134 = (-0.987289396019441,  0,                 0.158932842759347)
%c214 = ( 0.493644698009720, -0.855017697839831, 0.158932842759347)
%c243 = ( 0,                  0,                -1)

% Edge lengths of the inner tetrahedron
eB  = sqrt( (v(2,1)-v(1,1))^2 + (v(2,2)-v(1,2))^2 + (v(2,3)-v(1,3))^2 ); % eB = 1.903763289538843 % long edge with angle B
eAC = sqrt( (v(3,1)-v(2,1))^2 + (v(3,2)-v(2,2))^2 + (v(3,3)-v(2,3))^2 ); % eAC= 1.010545104217186 % short edge with angles A, C
     
% Angles of the inner tetrahedron faces
% edge formula would me much simpler :)
anA = acos( ( (v(2,1)-v(1,1))*(v(3,1)-v(1,1)) + (v(2,2)-v(1,2))*(v(3,2)-v(1,2)) + (v(2,3)-v(1,3))*(v(3,3)-v(1,3)) ) / ( sqrt( (v(2,1)-v(1,1))^2 + (v(2,2)-v(1,2))^2 + (v(2,3)-v(1,3))^2 )*sqrt( (v(3,1)-v(1,1))^2 + (v(3,2)-v(1,2))^2 + (v(3,3)-v(1,3))^2 ) ) );
anB = acos( ( (v(1,1)-v(2,1))*(v(3,1)-v(2,1)) + (v(1,2)-v(2,2))*(v(3,2)-v(2,2)) + (v(1,3)-v(2,3))*(v(3,3)-v(2,3)) ) / ( sqrt( (v(1,1)-v(2,1))^2 + (v(1,2)-v(2,2))^2 + (v(1,3)-v(2,3))^2 )*sqrt( (v(3,1)-v(2,1))^2 + (v(3,2)-v(2,2))^2 + (v(3,3)-v(2,3))^2 ) ) );
anC = acos( ( (v(4,1)-v(2,1))*(v(3,1)-v(2,1)) + (v(4,2)-v(2,2))*(v(3,2)-v(2,2)) + (v(4,3)-v(2,3))*(v(3,3)-v(2,3)) ) / ( sqrt( (v(4,1)-v(2,1))^2 + (v(4,2)-v(2,2))^2 + (v(4,3)-v(2,3))^2 )*sqrt( (v(3,1)-v(2,1))^2 + (v(3,2)-v(2,2))^2 + (v(3,3)-v(2,3))^2 ) ) ); 
%anA = 0.537252563377969; anAdeg = 30.782304414142388
%anB = 1.302170045105912; anBdeg = 74.608847792928799
%anC = pi/3;              anCdeg = 60

trg123 = eAC*sqrt(4*eB^2 - eAC^2)/4; % isosceles triangle (0.927421446334613) %ACAD OK
trg234 = (eAC^2)*sqrt(3)/4;          % equilateral triangle (0.442193180705835) %ACAD OK
%trg123a = eB*eB*sin(anA)/2          % 2 edges & angle formula
%trg123b = eAC*eB*sin(anB)/2         % 2 edges & angle formula
%trg234a = eAC*eAC*sin(anC)/2        % 2 edges & angle formula

trg123Planck = trg123/4; %eAC*sqrt(4*eB^2 - eAC^2)/16
trg234Planck = trg234/4; %sqrt(3)*(eAC^2)/16

NBHP = 3*trg123Planck+trg234Planck   % informational capacity of the inner tetrahedron
NBH  = pi;                           % informational capacity of the pi-bit black hole  
rtN  = NBHP/NBH %0.256594176525814

Afaces = 3*trg123 + trg234

% Laplacian of the inner tetrahedron
omit_1i = cot(anB);
omit_ij = (cot(anA)+cot(anC))/2;

L  = [ 3*omit_1i -omit_1i           -omit_1i           -omit_1i;
      -omit_1i    omit_1i+2*omit_ij -omit_ij           -omit_ij;
      -omit_1i   -omit_ij            omit_1i+2*omit_ij -omit_ij; 
      -omit_1i   -omit_ij           -omit_ij            omit_1i+2*omit_ij];
%eig(L) = [0 1.101119034602176 3.659346424371417 (x2)]
% The magnitude of the second smallest eigenvalue is related to the “mixing” properties of the graph (Chung, 1997).
% Namely, a smaller second eigenvalue indicates that a random walk on the graph will converge more rapidly to its stationary distribution.

% harmonic index
H=4*sum(eig(L)) %= 33.679247533380050
% the same as for faces   
%Hfaces1 = 12*(2*eB^2 + eAC^2)/eAC/sqrt(4*(eB^2) - (eAC^2)) + 12/sqrt(3)
%Hfaces2 = 3*(2*(eB^2) + (eAC^2))/trg123 + 3*(eAC^2)/trg234

%Voronoi tet
% 3*1.0771 + 1.2662 = 4.4975
      
%Ax = 1.868464859836944
%Ay = 1.868464859836944
%Az = 2.711985319745463      
%Aavg = (Ax+Ay+Az)/3 = 2.149638346473117
%Aavg/Afaces = 2/3
% -----------------------------------------------------------

% 4pi-bit black hole
% Spherical triangles angles
A = pi/3;              % = 60                 deg 
B = 1.476506070292930  % = 84.597566253231321 deg
C = 1.617889032100266  % = 92.698213259851002 deg
D = 0.188580513003934  % = 10.804867493537355 deg
E = 1.523703621489527  % = 87.301786740149012 deg
%E+C=pi, 2*B+D=pi, B+C=2*pi/3 + 1, 2*cos(C)+cos(1-C+2*pi/3)=0

% Cosine rule for angles
cosa = (cos(A)+cos(B)*cos(C))/(sin(B)*sin(C))
cosb = (cos(B)+cos(A)*cos(C))/(sin(A)*sin(C))
cosc = (cos(C)+cos(A)*cos(B))/(sin(A)*sin(B)) %=0
cosd = (cos(D)+cos(E)*cos(E))/(sin(E)*sin(E))
cose = (cos(E)+cos(D)*cos(E))/(sin(D)*sin(E))

% Spherical edges (angles in radians)
a = acos(cosa); % a = 1.049123144231259, adeg = 60.110328353945889 
b = acos(cosb); % b = 1.489078097122526, bdeg = 85.317890330492446
c = acos(cosc); % c = pi/2,              cdeg = 90
d = acos(cosd); % d = 0.163436459344741, ddeg = 9.364219339015140
e = acos(cose); % e = 1.049123144231253, edeg = 60.110328353945562
%a=e

% Cartesian coordinates on unit sphere (r=1)
p1 = [ 0                   0                   1     ]
p2 = [ 1                   0                   0     ]
p3 = [ cos(pi/3)*sin(b)    sin(pi/3)*sin(b)    cos(b)]
p4 = [ cos(2*pi/3)         sin(2*pi/3)         0     ]
p5 = [ cos(pi)*sin(b)      sin(pi)*sin(b)      cos(b)]
p6 = [ cos(4*pi/3)         sin(4*pi/3)         0     ]
p7 = [ cos(5*pi/3)*sin(b)  sin(5*pi/3)*sin(b)  cos(b)]
p8 = [ cos(pi/3)*sin(b)    sin(pi/3)*sin(b)   -cos(b)]
p9 = [ cos(pi)*sin(b)      sin(pi)*sin(b)     -cos(b)]
p10= [ cos(5*pi/3)*sin(b)  sin(5*pi/3)*sin(b) -cos(b)]
p11= [ 0                   0                  -1     ]
v=[p1;p2;p3;p4;p5;p6;p7;p8;p9;p10;p11];

%p1 =                   0,                  0,  1
%p2 =                   1,                  0,  0
%p3 =   0.498331461568344,  0.863135410446430,  0.081627309429082
%p4 =  -0.500000000000000,  0.866025403784439,                  0
%p5 =  -0.996662923136689,  0.000000000000000,  0.081627309429082
%p6 =  -0.500000000000000, -0.866025403784438,                  0
%p7 =   0.498331461568344, -0.863135410446430,  0.081627309429082
%p8 =   0.498331461568344,  0.863135410446430, -0.081627309429082
%p9 =  -0.996662923136689,  0.000000000000000, -0.081627309429082
%p10 =  0.498331461568344, -0.863135410446430, -0.081627309429082
%p11 =                  0,                  0,                 -1

% Voronoi vertices
c123   = cross(p2-p1, p3-p1)/norm(cross(p2-p1, p3-p1))
c134   = cross(p3-p1, p4-p1)/norm(cross(p3-p1, p4-p1))
c145   = cross(p4-p1, p5-p1)/norm(cross(p4-p1, p5-p1))
c156   = cross(p5-p1, p6-p1)/norm(cross(p5-p1, p6-p1))
c167   = cross(p6-p1, p7-p1)/norm(cross(p6-p1, p7-p1))
c172   = cross(p7-p1, p2-p1)/norm(cross(p7-p1, p2-p1))

c283   = cross(p8-p2, p3-p2 )/norm(cross(p8-p2, p3-p2 ))
c384   = cross(p8-p3, p4-p3 )/norm(cross(p8-p3, p4-p3 ))
c495   = cross(p9-p4, p5-p4 )/norm(cross(p9-p4, p5-p4 ))
c659   = cross(p5-p6, p9-p6 )/norm(cross(p5-p6, p9-p6 ))
c7610  = cross(p6-p7, p10-p7)/norm(cross(p6-p7, p10-p7))
c2710  = cross(p7-p2, p10-p2)/norm(cross(p7-p2, p10-p2))

c1182  = cross(p8-p11,  p2-p11 )/norm(cross(p8-p11,  p2-p11 ))
c1148  = cross(p4-p11,  p8-p11 )/norm(cross(p4-p11,  p8-p11 ))
c1194  = cross(p9-p11,  p4-p11 )/norm(cross(p9-p11,  p4-p11 ))
c1169  = cross(p6-p11,  p9-p11 )/norm(cross(p6-p11,  p9-p11 ))
c11106 = cross(p10-p11, p6-p11 )/norm(cross(p10-p11, p6-p11 ))
c11210 = cross(p2-p11,  p10-p11)/norm(cross(p2-p11,  p10-p11))

%c123  =   0.668627260652119,  0.325384653334613,  0.668627260652119
%c134  =  -0.052522254536692,  0.741740520054841,  0.668627260652119
%c145  =  -0.616105006115427,  0.416355866720229,  0.668627260652119
%c156  =  -0.616105006115427, -0.416355866720228,  0.668627260652119
%c167  =  -0.052522254536692, -0.741740520054841,  0.668627260652119
%c172  =   0.668627260652119, -0.325384653334613,  0.668627260652119
%c283  =   0.864574369493072,  0.502504885165963,                  0
%c384  =   0.002894811332970,  0.999995810024895,                  0
%c495  =  -0.867469180826042,  0.497490924858932,                  0
%c659  =  -0.867469180826042, -0.497490924858932,                  0
%c7610 =   0.002894811332970, -0.999995810024895,                  0
%c2710 =   0.864574369493072, -0.502504885165963,                  0
%c1182 =   0.668627260652119,  0.325384653334613, -0.668627260652119
%c1148 =  -0.052522254536692,  0.741740520054841, -0.668627260652119
%c1194 =  -0.616105006115427,  0.416355866720229, -0.668627260652119
%c1169 =  -0.616105006115427, -0.416355866720228, -0.668627260652119
%c11106=  -0.052522254536692, -0.741740520054841, -0.668627260652119
%c11210=   0.668627260652119, -0.325384653334613, -0.668627260652119

% Edge lengths of the inner polyhedron
eAE = sqrt( (p3(1)-p2(1))^2 + (p3(2)-p2(2))^2 + (p3(3)-p2(3))^2 ) % edge with angle A or E %ACAD OK
eB  = sqrt( (p3(1)-p1(1))^2 + (p3(2)-p1(2))^2 + (p3(3)-p1(3))^2 ) % edge with angle B %ACAD OK
eC  = sqrt( (p2(1)-p1(1))^2 + (p2(2)-p1(2))^2 + (p2(3)-p1(3))^2 ) % edge with angle C %ACAD OK
eD  = sqrt( (p8(1)-p3(1))^2 + (p8(2)-p3(2))^2 + (p8(3)-p3(3))^2 ) % edge with angle D %ACAD OK

% Angles of the inner polyhedron
%https://en.wikipedia.org/wiki/Triangle
anA = acos( (eB^2+eC^2-eAE^2)/(2*eB*eC) )
anB = acos( (eC^2+eAE^2-eB^2)/(2*eC*eAE) )
anC = acos( (eB^2+eAE^2-eC^2)/(2*eB*eAE) )
anD = acos( (2*eAE^2-eD^2)/(2*eAE^2) )
anE = acos( (eD^2)/(2*eAE*eD) )

trg123 = sqrt( (eC^2 + eAE^2 + eB^2)^2 - 2*(eC^4 + eAE^4 + eB^4) )/4       % 0.645453349901263 %ACAD OK
trg238 = eD*sqrt(4*eAE^2 - eD^2)/4 % 3 edges formula (isosceles triangle)  % 0.081491452573385 %ACAD OK

trg123Planck = trg123 % since r = lP
trg238Planck = trg238  

NBHP = 12*trg123Planck+6*trg238Planck  % informational capacity of the inner polyhedron
NBH  = 4*pi;                           % informational capacity of the 4pi-bit black hole  
rtN  = NBHP/NBH % 0.655271849522431

Afaces = 12*trg123 + 6*trg238 % 8.234388914255462

% cotan Laplacian
BB = cot(anB);
CC = cot(anC);
DD = cot(anD);
AE = (cot(anA)+cot(anE))/2;
d1 = 3*BB + 3*CC;
d2 = 2*CC + 4*AE;
d3 = BB + DD + 2*AE;
d4 = 2*CC + 4*AE;
d5 = BB + DD + 2*AE;
d6 = 2*CC + 4*AE;
d7 = BB + DD + 2*AE;
d8 = BB + DD + 2*AE;
d9 = BB + DD + 2*AE;
d10= BB + DD + 2*AE;
d11= 3*BB + 3*CC;

%     1   2   3   4   5   6   7   8   9   10   11
L = [ d1 -CC -BB -CC -BB -CC -BB  0   0   0    0  ; % 1
     -CC  d2 -AE  0   0   0  -AE -AE  0  -AE  -CC ; % 2  
     -BB -AE  d3 -AE  0   0   0  -DD  0   0    0  ; % 3
     -CC  0  -AE  d4 -AE  0   0  -AE -AE  0   -CC ; % 4
     -BB  0   0  -AE  d5 -AE  0   0  -DD  0    0  ; % 5
     -CC  0   0   0  -AE  d6 -AE  0  -AE -AE  -CC ; % 6
     -BB -AE  0   0   0  -AE  d7  0   0  -DD   0  ; % 7
     0   -AE -DD -AE  0   0   0   d8  0   0   -BB ; % 8
     0    0   0  -AE -DD -AE  0   0   d9  0   -BB ; % 9
     0   -AE  0   0   0  -AE -DD  0   0   d10 -BB ; % 10
     0   -CC  0  -CC  0  -CC  0  -BB -BB -BB   d11];% 11
%eig(L) = [0, 1.238424367281209 (x2), 2.278447744320292, 2.812918683486618, 3.401831584445386 (x2), 3.401831584445394, 4.159040427771492, 13.779710886934714 (x2), 13.832966302145932]

% harmonic index
H=4*sum(eig(L)) %=158.0065691210111
% the same as for faces   
Hfaces1 = 12*((eAE^2) + (eB^2)+ (eC^2))/trg123 + 6*(2*(eAE^2)+(eD^2))/trg238
% -----------------------------------------------------------
