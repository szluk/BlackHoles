% Draw plots for black hole informationless absortption and emission
% https://arxiv.org/abs/1910.11081 (Fig. 3)
% (c) Szymon Lukaszyk
% email: szymon@patent.pl
% licensed under MIT License.
% History
% v1:0 21.03.2021

clear all

RN_st  = 3.2; %  starting BH radius

% lP=1.616255*10^-35; %[m] Planck length
ali   = 4*pi^3+pi^2+pi
ali2  = -4*pi^3 - pi^2 - 2*pi;

% RN = R/lP % we make it dimensionless
RN_min = .1; % 1/2 is the minimum BH radius
RN_max = 8;

RN_dlt = (RN_max-RN_min)/100000;
RN     = RN_min:RN_dlt:RN_max;
RN2 = RN.^2; % = R^2/lP^2

nol = 100;
NAp_tab = zeros(nol, 3);
NAB_tab = zeros(nol, 3);

% NAP lines
RN2_st = RN_st^2;
for i=1:nol
  NAPl  = 64*pi^3/RN2_st + 32*pi^2 + 4*pi*RN2_st;
  NAl    = 4*pi*RN2_st;

  NAP_tab(i, 1) = sqrt(RN2_st);
  NAP_tab(i, 2) = NAPl;  
  NAP_tab(i, 3) = NAl;    
  RN2_st    = (NAPl)/(4*pi);  
end

% NAB lines
RN2_st = RN_st^2;
for i=1:nol
  NABl  = 64*pi^3/RN2_st - 32*pi^2 + 4*pi*RN2_st;
  NAl    = 4*pi*RN2_st;

  NAB_tab(i, 1) = sqrt(RN2_st);
  NAB_tab(i, 2) = NABl;  
  NAB_tab(i, 3) = NAl;    
  RN2_st    = (NABl)/(4*pi);  
end

NA  = 4*pi*RN2; % informational capacity 
NAP = (64*pi^3)./RN2 + 32*pi^2 + 4*pi*RN2; % after absorption
NAB = (64*pi^3)./RN2 - 32*pi^2 + 4*pi*RN2; % after emission

NAPN = floor(NAP);
NABN = floor(NAB);
NAN  = floor(NA);

figure
hold on
plot(RN, NAP, 'r')
plot(RN, NA,  'g')
plot(RN, NAB, 'b')

% NAp lines
for i=1:nol-1
  line([NAP_tab(i, 1)  NAP_tab(i, 1)  ], [NAP_tab(i, 3)  NAP_tab(i, 2)], 'Color', 'r')
  line([NAP_tab(i, 1)  NAP_tab(i+1, 1)], [NAP_tab(i, 2)  NAP_tab(i, 2)], 'Color', 'r')
end

% NAB lines
for i=1:nol-1
  line([NAB_tab(i, 1)  NAB_tab(i, 1)  ], [NAB_tab(i, 3)  NAB_tab(i, 2)])
  line([NAB_tab(i, 1)  NAB_tab(i+1, 1)], [NAB_tab(i, 2)  NAB_tab(i, 2)])
end

%plot(RN, NAPN, 'c')
%plot(RN, NABN, 'c')
%plot(RN, NAN,  'c')

axis([RN_min RN_max 0 2000])
