
function [t, STATES, ALGEBRAIC, CONSTANTS] = RasmussonUnparam(holding_p, holding_t, P1, P1_t, P2, P2_t)
    % This is the "main function".  In Matlab, things work best if you rename this function to match the filename.
   [t, STATES, ALGEBRAIC, CONSTANTS] = solveModel(holding_p, holding_t, P1, P1_t, P2, P2_t);
end

function [t, STATES, ALGEBRAIC, CONSTANTS] = solveModel(holding_p, holding_t, P1, P1_t, P2, P2_t)
    % Create ALGEBRAIC of correct size
    global algebraicVariableCount;  algebraicVariableCount = getAlgebraicVariableCount();
    
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts;

    % Set timespan to solve over 
    tspan = [0, P2_t];

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 1);

    % Solve model with ODE solver
    [t, STATES] = ode15s(@(t, STATES)computeRates(t, STATES, CONSTANTS, holding_p, holding_t, P1, P1_t, P2), tspan, INIT_STATES, options);

    % Compute algebraic variables
    [RATES, ALGEBRAIC] = computeRates(t, STATES, CONSTANTS, holding_p, holding_t, P1, P1_t, P2);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, t, holding_p, holding_t, P1, P1_t, P2);
end

function [RATES, ALGEBRAIC] = computeRates(t, STATES, CONSTANTS, holding_p, holding_t, P1, P1_t, P2)
    global algebraicVariableCount;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        ALGEBRAIC = zeros(1, algebraicVariableCount);
    else
        statesRowCount = statesSize(1);
        ALGEBRAIC = zeros(statesRowCount, algebraicVariableCount);
        RATES = zeros(statesRowCount, statesColumnCount);
    end
    
    % externally applied voltage (voltage clamp)
    ALGEBRAIC(:,72) = arrayfun(@(t) volt_clamp(t, holding_p, holding_t, P1, P1_t, P2), t);
    
    % A94; iKss if CONSTANTS(:,72)=0  or A110; sigma if CONSTANTS(:,72)=A110
    RATES(:,37) = CONSTANTS(:,72); 
    % A16; LTRPNCa
    RATES(:,7) =  CONSTANTS(:,24).*STATES(:,2).*(CONSTANTS(:,32) - STATES(:,7)) -  CONSTANTS(:,25).*STATES(:,7);
    % A17; HTEPNCa
    RATES(:,8) =  CONSTANTS(:,22).*STATES(:,2).*(CONSTANTS(:,33) - STATES(:,8)) -  CONSTANTS(:,23).*STATES(:,8);
    % A20; PO2
    RATES(:,10) =  CONSTANTS(:,37).*power(STATES(:,3), CONSTANTS(:,41)).*STATES(:,9) -  CONSTANTS(:,38).*STATES(:,10);
    % A21; PC2
    RATES(:,11) =  CONSTANTS(:,39).*STATES(:,9) -  CONSTANTS(:,40).*STATES(:,11);
    % A19; PC1
    ALGEBRAIC(:,2) = 1.00000 - (STATES(:,11)+STATES(:,9)+STATES(:,10));
    % A18; PO1
    RATES(:,9) = ( CONSTANTS(:,35).*power(STATES(:,3), CONSTANTS(:,42)).*ALGEBRAIC(:,2)+ CONSTANTS(:,38).*STATES(:,10)+ CONSTANTS(:,40).*STATES(:,11)) - ( CONSTANTS(:,36).*STATES(:,9)+ CONSTANTS(:,37).*power(STATES(:,3), CONSTANTS(:,41)).*STATES(:,9)+ CONSTANTS(:,39).*STATES(:,9));
    % A71; alpha_a
    ALGEBRAIC(:,5) =  0.180640.*exp( 0.0357700.*(ALGEBRAIC(:,72)+30.0000));
    % A72; beta_a
    ALGEBRAIC(:,15) =  0.395600.*exp(  - 0.0623700.*(ALGEBRAIC(:,72)+30.0000));
    % A69; ato_f
    RATES(:,29) =  ALGEBRAIC(:,5).*(1.00000 - STATES(:,29)) -  ALGEBRAIC(:,15).*STATES(:,29);
    % A73; alpha_i
    ALGEBRAIC(:,6) = ( 0.000152000.*exp( - (ALGEBRAIC(:,72)+13.5000)./7.00000))./( 0.00670830.*exp( - (ALGEBRAIC(:,72)+33.5000)./7.00000)+1.00000);
    % A74; beta_i
    ALGEBRAIC(:,16) = ( 0.000950000.*exp((ALGEBRAIC(:,72)+33.5000)./7.00000))./( 0.0513350.*exp((ALGEBRAIC(:,72)+33.5000)./7.00000)+1.00000);
    % A70; ito_f
    RATES(:,30) =  ALGEBRAIC(:,6).*(1.00000 - STATES(:,30)) -  ALGEBRAIC(:,16).*STATES(:,30);
    % A78; a_ss
    ALGEBRAIC(:,7) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,72)+22.5000)./7.70000));
    % A80; tau_tas
    ALGEBRAIC(:,17) =  0.493000.*exp(  - 0.0629000.*ALGEBRAIC(:,72))+2.05800;
    % A76; ato_s
    RATES(:,31) = (ALGEBRAIC(:,7) - STATES(:,31))./ALGEBRAIC(:,17);
    % A79; i_ss
    ALGEBRAIC(:,8) = 1.00000./(1.00000+exp((ALGEBRAIC(:,72)+45.2000)./5.70000));
    % A81; tau_tis
    ALGEBRAIC(:,18) = 270.000+1050.00./(1.00000+exp((ALGEBRAIC(:,72)+45.2000)./5.70000));
    % A77; ito_s
    RATES(:,32) = (ALGEBRAIC(:,8) - STATES(:,32))./ALGEBRAIC(:,18);
    % A85; alpha_n
    ALGEBRAIC(:,9) = ( 4.81333e-06.*(ALGEBRAIC(:,72)+26.5000))./(1.00000 - exp(  - 0.128000.*(ALGEBRAIC(:,72)+26.5000)));
    % A86; beta_n
    ALGEBRAIC(:,19) =  9.53333e-05.*exp(  - 0.0380000.*(ALGEBRAIC(:,72)+26.5000));
    % A84; nKs
    RATES(:,33) =  ALGEBRAIC(:,9).*(1.00000 - STATES(:,33)) -  ALGEBRAIC(:,19).*STATES(:,33);
    % A90; tau_aur
    ALGEBRAIC(:,20) =  0.493000.*exp(  - 0.0629000.*ALGEBRAIC(:,72))+2.05800;
    % A88; aur
    RATES(:,34) = (ALGEBRAIC(:,7) - STATES(:,34))./ALGEBRAIC(:,20);
    % A91; tau_iur
    ALGEBRAIC(:,21) = 1200.00 - 170.000./(1.00000+exp((ALGEBRAIC(:,72)+45.2000)./5.70000));
    % A89; iur
    RATES(:,35) = (ALGEBRAIC(:,8) - STATES(:,35))./ALGEBRAIC(:,21);
    % A95; tau_Kss
    ALGEBRAIC(:,22) =  39.3000.*exp(  - 0.0862000.*ALGEBRAIC(:,72))+13.1700;
    % A93; aKss
    RATES(:,36) = (ALGEBRAIC(:,7) - STATES(:,36))./ALGEBRAIC(:,22);
    % A104; alpha_a1
    ALGEBRAIC(:,11) =  0.0137330.*exp( 0.0381980.*ALGEBRAIC(:,72));
    % A105; beta_a1
    ALGEBRAIC(:,24) =  6.89000e-05.*exp(  - 0.0417800.*ALGEBRAIC(:,72));
    % A99; CK2
    RATES(:,40) = ( CONSTANTS(:,65).*STATES(:,39)+ ALGEBRAIC(:,24).*STATES(:,38)) - ( CONSTANTS(:,64).*STATES(:,40)+ ALGEBRAIC(:,11).*STATES(:,40));
    % A24; C1
    ALGEBRAIC(:,3) = 1.00000 - (STATES(:,12)+STATES(:,13)+STATES(:,14)+STATES(:,15)+STATES(:,16)+STATES(:,17)+STATES(:,18));
    % A31; alpha
    ALGEBRAIC(:,13) = ( 0.400000.*exp((ALGEBRAIC(:,72)+12.0000)./10.0000).*((1.00000+ 0.700000.*exp( - power(ALGEBRAIC(:,72)+40.0000, 2.00000)./10.0000)) -  0.750000.*exp( - power(ALGEBRAIC(:,72)+20.0000, 2.00000)./400.000)))./(1.00000+ 0.120000.*exp((ALGEBRAIC(:,72)+12.0000)./10.0000));
    % A32; beta
    ALGEBRAIC(:,26) =  0.0500000.*exp( - (ALGEBRAIC(:,72)+12.0000)./13.0000);
    % A25; C2
    RATES(:,13) = ( 4.00000.*ALGEBRAIC(:,13).*ALGEBRAIC(:,3)+ 2.00000.*ALGEBRAIC(:,26).*STATES(:,14)) - ( ALGEBRAIC(:,26).*STATES(:,13)+ 3.00000.*ALGEBRAIC(:,13).*STATES(:,13));
    % A26; C3
    RATES(:,14) = ( 3.00000.*ALGEBRAIC(:,13).*STATES(:,13)+ 3.00000.*ALGEBRAIC(:,26).*STATES(:,15)) - ( 2.00000.*ALGEBRAIC(:,26).*STATES(:,14)+ 2.00000.*ALGEBRAIC(:,13).*STATES(:,14));
    % A97; CK0
    ALGEBRAIC(:,10) = 1.00000 - (STATES(:,39)+STATES(:,40)+STATES(:,38)+STATES(:,41));
    % A102; alpha_a0
    ALGEBRAIC(:,23) =  0.0223480.*exp( 0.0117600.*ALGEBRAIC(:,72));
    % A103; beta_a0
    ALGEBRAIC(:,28) =  0.0470020.*exp(  - 0.0631000.*ALGEBRAIC(:,72));
    % A98; CK1
    RATES(:,39) = ( ALGEBRAIC(:,23).*ALGEBRAIC(:,10)+ CONSTANTS(:,64).*STATES(:,40)) - ( ALGEBRAIC(:,28).*STATES(:,39)+ CONSTANTS(:,65).*STATES(:,39));
    % A106; alpha_i
    ALGEBRAIC(:,29) =  0.0908210.*exp( 0.0233910.*(ALGEBRAIC(:,72)+5.00000));
    % A107; beta_i
    ALGEBRAIC(:,33) =  0.00649700.*exp(  - 0.0326800.*(ALGEBRAIC(:,72)+5.00000));
    % A100; OK
    RATES(:,38) = ( ALGEBRAIC(:,11).*STATES(:,40)+ ALGEBRAIC(:,33).*STATES(:,41)) - ( ALGEBRAIC(:,24).*STATES(:,38)+ ALGEBRAIC(:,29).*STATES(:,38));
    % A101; IK
    RATES(:,41) =  ALGEBRAIC(:,29).*STATES(:,38) -  ALGEBRAIC(:,33).*STATES(:,41);
    % A33; gamma
    ALGEBRAIC(:,31) = ( CONSTANTS(:,46).*STATES(:,3))./(CONSTANTS(:,47)+STATES(:,3));
    % A34; K_pcf
    ALGEBRAIC(:,35) =  13.0000.*(1.00000 - exp( - power(ALGEBRAIC(:,72)+14.5000, 2.00000)./100.000));
    % A23; O
    RATES(:,12) = ( ALGEBRAIC(:,13).*STATES(:,15)+ CONSTANTS(:,45).*STATES(:,16)+ 0.00100000.*( ALGEBRAIC(:,13).*STATES(:,17) -  ALGEBRAIC(:,35).*STATES(:,12))) - ( 4.00000.*ALGEBRAIC(:,26).*STATES(:,12)+ ALGEBRAIC(:,31).*STATES(:,12));
    % A27; C4
    RATES(:,15) = ( 2.00000.*ALGEBRAIC(:,13).*STATES(:,14)+ 4.00000.*ALGEBRAIC(:,26).*STATES(:,12)+ 0.0100000.*( 4.00000.*CONSTANTS(:,45).*ALGEBRAIC(:,26).*STATES(:,16) -  ALGEBRAIC(:,13).*ALGEBRAIC(:,31).*STATES(:,15))+ 0.00200000.*( 4.00000.*ALGEBRAIC(:,26).*STATES(:,17) -  ALGEBRAIC(:,35).*STATES(:,15))+ 4.00000.*ALGEBRAIC(:,26).*CONSTANTS(:,45).*STATES(:,18)) - ( 3.00000.*ALGEBRAIC(:,26).*STATES(:,15)+ ALGEBRAIC(:,13).*STATES(:,15)+ 1.00000.*ALGEBRAIC(:,31).*ALGEBRAIC(:,35).*STATES(:,15));
    % A28; I1
    RATES(:,16) = ( ALGEBRAIC(:,31).*STATES(:,12)+ 0.00100000.*( ALGEBRAIC(:,13).*STATES(:,18) -  ALGEBRAIC(:,35).*STATES(:,16))+ 0.0100000.*( ALGEBRAIC(:,13).*ALGEBRAIC(:,31).*STATES(:,15) -  4.00000.*ALGEBRAIC(:,26).*ALGEBRAIC(:,35).*STATES(:,16))) -  CONSTANTS(:,45).*STATES(:,16);
    % A29; I2
    RATES(:,17) = ( 0.00100000.*( ALGEBRAIC(:,35).*STATES(:,12) -  ALGEBRAIC(:,13).*STATES(:,17))+ CONSTANTS(:,45).*STATES(:,18)+ 0.00200000.*( ALGEBRAIC(:,35).*STATES(:,15) -  4.00000.*ALGEBRAIC(:,26).*STATES(:,17))) -  ALGEBRAIC(:,31).*STATES(:,17);
    % A30; I3
    RATES(:,18) = ( 0.00100000.*( ALGEBRAIC(:,35).*STATES(:,16) -  ALGEBRAIC(:,13).*STATES(:,18))+ ALGEBRAIC(:,31).*STATES(:,17)+ 1.00000.*ALGEBRAIC(:,31).*ALGEBRAIC(:,35).*STATES(:,15)) - ( 4.00000.*ALGEBRAIC(:,26).*CONSTANTS(:,45).*STATES(:,18)+ CONSTANTS(:,45).*STATES(:,18));
    % A8; B_JSR
    ALGEBRAIC(:,30) = power(1.00000+( CONSTANTS(:,19).*CONSTANTS(:,21))./power(CONSTANTS(:,21)+STATES(:,4), 2.00000),  - 1.00000);
    % A9; J_rel
    ALGEBRAIC(:,34) =  CONSTANTS(:,26).*(STATES(:,9)+STATES(:,10)).*(STATES(:,4) - STATES(:,3)).*STATES(:,6);
    % A10; J_tr
    ALGEBRAIC(:,37) = (STATES(:,5) - STATES(:,4))./CONSTANTS(:,27);
    % A4; CaJSR
    RATES(:,4) =  ALGEBRAIC(:,30).*(ALGEBRAIC(:,37) - ALGEBRAIC(:,34));
    % A12; J_leak
    ALGEBRAIC(:,41) =  CONSTANTS(:,28).*(STATES(:,5) - STATES(:,2));
    % A13; J_up
    ALGEBRAIC(:,43) = ( CONSTANTS(:,30).*power(STATES(:,2), 2.00000))./(power(CONSTANTS(:,31), 2.00000)+power(STATES(:,2), 2.00000));
    % A5; CaNSR
    RATES(:,5) = ( (ALGEBRAIC(:,43) - ALGEBRAIC(:,41)).*CONSTANTS(:,2))./CONSTANTS(:,4) - ( ALGEBRAIC(:,37).*CONSTANTS(:,3))./CONSTANTS(:,4);
    % A42; CNa3
    ALGEBRAIC(:,4) = 1.00000 - (STATES(:,20)+STATES(:,21)+STATES(:,22)+STATES(:,25)+STATES(:,23)+STATES(:,24)+STATES(:,26)+STATES(:,27));
    % A51; alpha_Na11
    ALGEBRAIC(:,14) = 3.80200./( 0.102700.*exp( - (ALGEBRAIC(:,72)+2.50000)./17.0000)+ 0.200000.*exp( - (ALGEBRAIC(:,72)+2.50000)./150.000));
    % A54; beta_Na11
    ALGEBRAIC(:,36) =  0.191700.*exp( - (ALGEBRAIC(:,72)+2.50000)./20.3000);
    % A52; alpha_Na12
    ALGEBRAIC(:,27) = 3.80200./( 0.102700.*exp( - (ALGEBRAIC(:,72)+2.50000)./15.0000)+ 0.230000.*exp( - (ALGEBRAIC(:,72)+2.50000)./150.000));
    % A55; beta_Na12
    ALGEBRAIC(:,38) =  0.200000.*exp( - (ALGEBRAIC(:,72) - 2.50000)./20.3000);
    % A57; alpha_Na3
    ALGEBRAIC(:,42) =  7.00000e-07.*exp( - (ALGEBRAIC(:,72)+7.00000)./7.70000);
    % A58; beta_Na3
    ALGEBRAIC(:,44) = 0.00840000+ 2.00000e-05.*(ALGEBRAIC(:,72)+7.00000);
    % A43; CNa2
    RATES(:,22) = ( ALGEBRAIC(:,14).*ALGEBRAIC(:,4)+ ALGEBRAIC(:,38).*STATES(:,21)+ ALGEBRAIC(:,42).*STATES(:,26)) - ( ALGEBRAIC(:,36).*STATES(:,22)+ ALGEBRAIC(:,27).*STATES(:,22)+ ALGEBRAIC(:,44).*STATES(:,22));
    % A53; alpha_Na13
    ALGEBRAIC(:,32) = 3.80200./( 0.102700.*exp( - (ALGEBRAIC(:,72)+2.50000)./12.0000)+ 0.250000.*exp( - (ALGEBRAIC(:,72)+2.50000)./150.000));
    % A56; beta_Na13
    ALGEBRAIC(:,40) =  0.220000.*exp( - (ALGEBRAIC(:,72) - 7.50000)./20.3000);
    % A44; CNa1
    RATES(:,21) = ( ALGEBRAIC(:,27).*STATES(:,22)+ ALGEBRAIC(:,40).*STATES(:,20)+ ALGEBRAIC(:,42).*STATES(:,25)) - ( ALGEBRAIC(:,38).*STATES(:,21)+ ALGEBRAIC(:,32).*STATES(:,21)+ ALGEBRAIC(:,44).*STATES(:,21));
    % A49; ICNa2
    RATES(:,26) = ( ALGEBRAIC(:,14).*STATES(:,27)+ ALGEBRAIC(:,38).*STATES(:,25)+ ALGEBRAIC(:,44).*STATES(:,22)) - ( ALGEBRAIC(:,36).*STATES(:,26)+ ALGEBRAIC(:,27).*STATES(:,26)+ ALGEBRAIC(:,42).*STATES(:,26));
    % A50; ICNa3
    RATES(:,27) = ( ALGEBRAIC(:,36).*STATES(:,26)+ ALGEBRAIC(:,44).*ALGEBRAIC(:,4)) - ( ALGEBRAIC(:,14).*STATES(:,27)+ ALGEBRAIC(:,42).*STATES(:,27));
    % A22; I_CaL
    ALGEBRAIC(:,47) =  CONSTANTS(:,44).*STATES(:,12).*(ALGEBRAIC(:,72) - CONSTANTS(:,43));
    % A7; B_ss
    ALGEBRAIC(:,25) = power(1.00000+( CONSTANTS(:,18).*CONSTANTS(:,20))./power(CONSTANTS(:,20)+STATES(:,3), 2.00000),  - 1.00000);
    % A11; J_xfer
    ALGEBRAIC(:,39) = (STATES(:,3) - STATES(:,2))./CONSTANTS(:,29);
    % A3; Cass
    RATES(:,3) =  ALGEBRAIC(:,25).*(( ALGEBRAIC(:,34).*CONSTANTS(:,3))./CONSTANTS(:,5) - (( ALGEBRAIC(:,39).*CONSTANTS(:,2))./CONSTANTS(:,5)+( ALGEBRAIC(:,47).*CONSTANTS(:,6).*CONSTANTS(:,1))./( 2.00000.*CONSTANTS(:,5).*CONSTANTS(:,12))));
    % A15; PRyR
    RATES(:,6) =   - 0.0400000.*STATES(:,6) -  (( 0.100000.*ALGEBRAIC(:,47))./CONSTANTS(:,34)).*exp( - power(ALGEBRAIC(:,72) - 5.00000, 2.00000)./648.000);
    % A59; alpha_Na2
    ALGEBRAIC(:,46) = 1.00000./( 0.188495.*exp( - (ALGEBRAIC(:,72)+7.00000)./16.6000)+0.393956);
    % A60; beta_Na2
    ALGEBRAIC(:,48) = ( ALGEBRAIC(:,32).*ALGEBRAIC(:,46).*ALGEBRAIC(:,42))./( ALGEBRAIC(:,40).*ALGEBRAIC(:,44));
    % A45; ONa
    RATES(:,20) = ( ALGEBRAIC(:,32).*STATES(:,21)+ ALGEBRAIC(:,48).*STATES(:,25)) - ( ALGEBRAIC(:,40).*STATES(:,20)+ ALGEBRAIC(:,46).*STATES(:,20));
    % A61; alpha_Na4
    ALGEBRAIC(:,50) = ALGEBRAIC(:,46)./1000.00;
    % A62; beta_Na4
    ALGEBRAIC(:,52) = ALGEBRAIC(:,42);
    % A46; IFNa
    RATES(:,25) = ( ALGEBRAIC(:,46).*STATES(:,20)+ ALGEBRAIC(:,44).*STATES(:,21)+ ALGEBRAIC(:,52).*STATES(:,23)+ ALGEBRAIC(:,27).*STATES(:,26)) - ( ALGEBRAIC(:,48).*STATES(:,25)+ ALGEBRAIC(:,42).*STATES(:,25)+ ALGEBRAIC(:,50).*STATES(:,25)+ ALGEBRAIC(:,38).*STATES(:,25));
    % A35; I_p(Ca)
    ALGEBRAIC(:,49) = ( CONSTANTS(:,48).*power(STATES(:,2), 2.00000))./(power(CONSTANTS(:,49), 2.00000)+power(STATES(:,2), 2.00000));
    % A36; I_NaCa
    ALGEBRAIC(:,51) =  (( (( (( CONSTANTS(:,50).*1.00000)./(power(CONSTANTS(:,51), 3.00000)+power(CONSTANTS(:,8), 3.00000))).*1.00000)./(CONSTANTS(:,52)+CONSTANTS(:,9))).*1.00000)./(1.00000+ CONSTANTS(:,53).*exp(( (CONSTANTS(:,54) - 1.00000).*ALGEBRAIC(:,72).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11))))).*( exp(( CONSTANTS(:,54).*ALGEBRAIC(:,72).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11))).*power(STATES(:,19), 3.00000).*CONSTANTS(:,9) -  exp(( (CONSTANTS(:,54) - 1.00000).*ALGEBRAIC(:,72).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11))).*power(CONSTANTS(:,8), 3.00000).*STATES(:,2));
    % A38; E_CaN
    ALGEBRAIC(:,53) =  (( CONSTANTS(:,10).*CONSTANTS(:,11))./( 2.00000.*CONSTANTS(:,12))).*log(CONSTANTS(:,9)./STATES(:,2));
    % A37; I_Cab
    ALGEBRAIC(:,55) =  CONSTANTS(:,55).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,53));
    % A6; Bi
    ALGEBRAIC(:,12) = power(1.00000+( CONSTANTS(:,18).*CONSTANTS(:,20))./power(CONSTANTS(:,20)+STATES(:,2), 2.00000),  - 1.00000);
    % A14; J_trpn
    ALGEBRAIC(:,45) = ( CONSTANTS(:,22).*STATES(:,2).*(CONSTANTS(:,33) - STATES(:,8))+ CONSTANTS(:,24).*STATES(:,2).*(CONSTANTS(:,32) - STATES(:,7))) - ( CONSTANTS(:,23).*STATES(:,8)+ CONSTANTS(:,25).*STATES(:,7));
    % A2; Cai
    RATES(:,2) =  ALGEBRAIC(:,12).*((ALGEBRAIC(:,41)+ALGEBRAIC(:,39)) - (ALGEBRAIC(:,43)+ALGEBRAIC(:,45)+( ((ALGEBRAIC(:,55)+ALGEBRAIC(:,49)) -  2.00000.*ALGEBRAIC(:,51)).*CONSTANTS(:,6).*CONSTANTS(:,1))./( 2.00000.*CONSTANTS(:,2).*CONSTANTS(:,12))));
    % A63; alpha_Na5
    ALGEBRAIC(:,54) = ALGEBRAIC(:,46)./95000.0;
    % A64; beta_Na5
    ALGEBRAIC(:,56) = ALGEBRAIC(:,42)./50.0000;
    % A47; I1Na
    RATES(:,23) = ( ALGEBRAIC(:,50).*STATES(:,25)+ ALGEBRAIC(:,56).*STATES(:,24)) - ( ALGEBRAIC(:,52).*STATES(:,23)+ ALGEBRAIC(:,54).*STATES(:,23));
    % A48; I2Na
    RATES(:,24) =  ALGEBRAIC(:,54).*STATES(:,23) -  ALGEBRAIC(:,56).*STATES(:,24);
    % A41; E_Na
    ALGEBRAIC(:,57) =  (( CONSTANTS(:,10).*CONSTANTS(:,11))./CONSTANTS(:,12)).*log(( 0.900000.*CONSTANTS(:,8)+ 0.100000.*CONSTANTS(:,7))./( 0.900000.*STATES(:,19)+ 0.100000.*STATES(:,28)));
    % A40; I_Na
    ALGEBRAIC(:,58) =  CONSTANTS(:,56).*STATES(:,20).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,57));
    % A65; I_Nab
    ALGEBRAIC(:,59) =  CONSTANTS(:,57).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,57));
    % A109; f_NaK
    ALGEBRAIC(:,68) = 1.00000./(1.00000+ 0.124500.*exp((  - 0.100000.*ALGEBRAIC(:,72).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11)))+ 0.0365000.*CONSTANTS(:,72).*exp((  - ALGEBRAIC(:,72).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11))));
    % A108; I_NaK
    ALGEBRAIC(:,69) = ( (( CONSTANTS(:,66).*ALGEBRAIC(:,68).*1.00000)./(1.00000+power(CONSTANTS(:,67)./STATES(:,19), 1.50000))).*CONSTANTS(:,7))./(CONSTANTS(:,7)+CONSTANTS(:,68));
    % A39; Nai
    RATES(:,19) = (  - (ALGEBRAIC(:,58)+ALGEBRAIC(:,59)+ 3.00000.*ALGEBRAIC(:,69)+ 3.00000.*ALGEBRAIC(:,51)).*CONSTANTS(:,6).*CONSTANTS(:,1))./( CONSTANTS(:,2).*CONSTANTS(:,12));
    % A68; E_K
    ALGEBRAIC(:,60) =  (( CONSTANTS(:,10).*CONSTANTS(:,11))./CONSTANTS(:,12)).*log(CONSTANTS(:,7)./STATES(:,28));
    % A67; I_Kto,f  
    ALGEBRAIC(:,61) =  CONSTANTS(:,58).*power(STATES(:,29), 3.00000).*STATES(:,30).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60));
    % A75; I_Kto,s
    ALGEBRAIC(:,62) =  CONSTANTS(:,59).*STATES(:,31).*STATES(:,32).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60));
    % A82; I_K1
    ALGEBRAIC(:,63) = ( (( 0.293800.*CONSTANTS(:,7))./(CONSTANTS(:,7)+210.000)).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60)))./(1.00000+exp( 0.0896000.*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60))));
    % A83; I_Ks
    ALGEBRAIC(:,64) =  CONSTANTS(:,60).*power(STATES(:,33), 2.00000).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60));
    % A87; I_kUR
    ALGEBRAIC(:,65) =  CONSTANTS(:,61).*STATES(:,34).*STATES(:,35).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60));
    % A92; I_Kss
    ALGEBRAIC(:,66) =  CONSTANTS(:,62).*STATES(:,36).*STATES(:,37).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60));
    % A96; I_Kr
    ALGEBRAIC(:,67) =  CONSTANTS(:,63).*STATES(:,38).*(ALGEBRAIC(:,72) -  (( CONSTANTS(:,10).*CONSTANTS(:,11))./CONSTANTS(:,12)).*log(( 0.980000.*CONSTANTS(:,7)+ 0.0200000.*CONSTANTS(:,8))./( 0.980000.*STATES(:,28)+ 0.0200000.*STATES(:,19))));
    % A66; Ki
    RATES(:,28) = (  - ((ALGEBRAIC(:,61)+ALGEBRAIC(:,62)+ALGEBRAIC(:,63)+ALGEBRAIC(:,64)+ALGEBRAIC(:,66)+ALGEBRAIC(:,65)+ALGEBRAIC(:,67)) -  2.00000.*ALGEBRAIC(:,69)).*CONSTANTS(:,6).*CONSTANTS(:,1))./( CONSTANTS(:,2).*CONSTANTS(:,12));
    % I_stim
    ALGEBRAIC(:,1) = piecewise({t>=CONSTANTS(:,13)&t<=CONSTANTS(:,14)&(t - CONSTANTS(:,13)) -  floor((t - CONSTANTS(:,13))./CONSTANTS(:,15)).*CONSTANTS(:,15)<=CONSTANTS(:,16), CONSTANTS(:,17) }, 0.00000);
    % A112; O_ClCa
    ALGEBRAIC(:,70) = 0.200000./(1.00000+exp( - (ALGEBRAIC(:,72) - 46.7000)./7.80000));
    % A111; i_ClCa
    ALGEBRAIC(:,71) =  (( CONSTANTS(:,69).*ALGEBRAIC(:,70).*STATES(:,2))./(STATES(:,2)+CONSTANTS(:,71))).*(ALGEBRAIC(:,72) - CONSTANTS(:,70));
    % A1; V
    RATES(:,1) =  0;
    RATES = RATES';
end

function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, t, holding_p, holding_t, P1, P1_t, P2)
    ALGEBRAIC(:,72) = arrayfun(@(t) volt_clamp(t, holding_p, holding_t, P1, P1_t, P2), t);

    ALGEBRAIC(:,2) = 1.00000 - (STATES(:,11)+STATES(:,9)+STATES(:,10));
    ALGEBRAIC(:,5) =  0.180640.*exp( 0.0357700.*(ALGEBRAIC(:,72)+30.0000));
    ALGEBRAIC(:,15) =  0.395600.*exp(  - 0.0623700.*(ALGEBRAIC(:,72)+30.0000));
    ALGEBRAIC(:,6) = ( 0.000152000.*exp( - (ALGEBRAIC(:,72)+13.5000)./7.00000))./( 0.00670830.*exp( - (ALGEBRAIC(:,72)+33.5000)./7.00000)+1.00000);
    ALGEBRAIC(:,16) = ( 0.000950000.*exp((ALGEBRAIC(:,72)+33.5000)./7.00000))./( 0.0513350.*exp((ALGEBRAIC(:,72)+33.5000)./7.00000)+1.00000);
    ALGEBRAIC(:,7) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,72)+22.5000)./7.70000));
    ALGEBRAIC(:,17) =  0.493000.*exp(  - 0.0629000.*ALGEBRAIC(:,72))+2.05800;
    ALGEBRAIC(:,8) = 1.00000./(1.00000+exp((ALGEBRAIC(:,72)+45.2000)./5.70000));
    ALGEBRAIC(:,18) = 270.000+1050.00./(1.00000+exp((ALGEBRAIC(:,72)+45.2000)./5.70000));
    ALGEBRAIC(:,9) = ( 4.81333e-06.*(ALGEBRAIC(:,72)+26.5000))./(1.00000 - exp(  - 0.128000.*(ALGEBRAIC(:,72)+26.5000)));
    ALGEBRAIC(:,19) =  9.53333e-05.*exp(  - 0.0380000.*(ALGEBRAIC(:,72)+26.5000));
    ALGEBRAIC(:,20) =  0.493000.*exp(  - 0.0629000.*ALGEBRAIC(:,72))+2.05800;
    ALGEBRAIC(:,21) = 1200.00 - 170.000./(1.00000+exp((ALGEBRAIC(:,72)+45.2000)./5.70000));
    ALGEBRAIC(:,22) =  39.3000.*exp(  - 0.0862000.*ALGEBRAIC(:,72))+13.1700;
    ALGEBRAIC(:,11) =  0.0137330.*exp( 0.0381980.*ALGEBRAIC(:,72));
    ALGEBRAIC(:,24) =  6.89000e-05.*exp(  - 0.0417800.*ALGEBRAIC(:,72));
    ALGEBRAIC(:,3) = 1.00000 - (STATES(:,12)+STATES(:,13)+STATES(:,14)+STATES(:,15)+STATES(:,16)+STATES(:,17)+STATES(:,18));
    ALGEBRAIC(:,13) = ( 0.400000.*exp((ALGEBRAIC(:,72)+12.0000)./10.0000).*((1.00000+ 0.700000.*exp( - power(ALGEBRAIC(:,72)+40.0000, 2.00000)./10.0000)) -  0.750000.*exp( - power(ALGEBRAIC(:,72)+20.0000, 2.00000)./400.000)))./(1.00000+ 0.120000.*exp((ALGEBRAIC(:,72)+12.0000)./10.0000));
    ALGEBRAIC(:,26) =  0.0500000.*exp( - (ALGEBRAIC(:,72)+12.0000)./13.0000);
    ALGEBRAIC(:,10) = 1.00000 - (STATES(:,39)+STATES(:,40)+STATES(:,38)+STATES(:,41));
    ALGEBRAIC(:,23) =  0.0223480.*exp( 0.0117600.*ALGEBRAIC(:,72));
    ALGEBRAIC(:,28) =  0.0470020.*exp(  - 0.0631000.*ALGEBRAIC(:,72));
    ALGEBRAIC(:,29) =  0.0908210.*exp( 0.0233910.*(ALGEBRAIC(:,72)+5.00000));
    ALGEBRAIC(:,33) =  0.00649700.*exp(  - 0.0326800.*(ALGEBRAIC(:,72)+5.00000));
    ALGEBRAIC(:,31) = ( CONSTANTS(:,46).*STATES(:,3))./(CONSTANTS(:,47)+STATES(:,3));
    ALGEBRAIC(:,35) =  13.0000.*(1.00000 - exp( - power(ALGEBRAIC(:,72)+14.5000, 2.00000)./100.000));
    ALGEBRAIC(:,30) = power(1.00000+( CONSTANTS(:,19).*CONSTANTS(:,21))./power(CONSTANTS(:,21)+STATES(:,4), 2.00000),  - 1.00000);
    ALGEBRAIC(:,34) =  CONSTANTS(:,26).*(STATES(:,9)+STATES(:,10)).*(STATES(:,4) - STATES(:,3)).*STATES(:,6);
    ALGEBRAIC(:,37) = (STATES(:,5) - STATES(:,4))./CONSTANTS(:,27);
    ALGEBRAIC(:,41) =  CONSTANTS(:,28).*(STATES(:,5) - STATES(:,2));
    ALGEBRAIC(:,43) = ( CONSTANTS(:,30).*power(STATES(:,2), 2.00000))./(power(CONSTANTS(:,31), 2.00000)+power(STATES(:,2), 2.00000));
    ALGEBRAIC(:,4) = 1.00000 - (STATES(:,20)+STATES(:,21)+STATES(:,22)+STATES(:,25)+STATES(:,23)+STATES(:,24)+STATES(:,26)+STATES(:,27));
    ALGEBRAIC(:,14) = 3.80200./( 0.102700.*exp( - (ALGEBRAIC(:,72)+2.50000)./17.0000)+ 0.200000.*exp( - (ALGEBRAIC(:,72)+2.50000)./150.000));
    ALGEBRAIC(:,36) =  0.191700.*exp( - (ALGEBRAIC(:,72)+2.50000)./20.3000);
    ALGEBRAIC(:,27) = 3.80200./( 0.102700.*exp( - (ALGEBRAIC(:,72)+2.50000)./15.0000)+ 0.230000.*exp( - (ALGEBRAIC(:,72)+2.50000)./150.000));
    ALGEBRAIC(:,38) =  0.200000.*exp( - (ALGEBRAIC(:,72) - 2.50000)./20.3000);
    ALGEBRAIC(:,42) =  7.00000e-07.*exp( - (ALGEBRAIC(:,72)+7.00000)./7.70000);
    ALGEBRAIC(:,44) = 0.00840000+ 2.00000e-05.*(ALGEBRAIC(:,72)+7.00000);
    ALGEBRAIC(:,32) = 3.80200./( 0.102700.*exp( - (ALGEBRAIC(:,72)+2.50000)./12.0000)+ 0.250000.*exp( - (ALGEBRAIC(:,72)+2.50000)./150.000));
    ALGEBRAIC(:,40) =  0.220000.*exp( - (ALGEBRAIC(:,72) - 7.50000)./20.3000);
    ALGEBRAIC(:,47) =  CONSTANTS(:,44).*STATES(:,12).*(ALGEBRAIC(:,72) - CONSTANTS(:,43));
    ALGEBRAIC(:,25) = power(1.00000+( CONSTANTS(:,18).*CONSTANTS(:,20))./power(CONSTANTS(:,20)+STATES(:,3), 2.00000),  - 1.00000);
    ALGEBRAIC(:,39) = (STATES(:,3) - STATES(:,2))./CONSTANTS(:,29);
    ALGEBRAIC(:,46) = 1.00000./( 0.188495.*exp( - (ALGEBRAIC(:,72)+7.00000)./16.6000)+0.393956);
    ALGEBRAIC(:,48) = ( ALGEBRAIC(:,32).*ALGEBRAIC(:,46).*ALGEBRAIC(:,42))./( ALGEBRAIC(:,40).*ALGEBRAIC(:,44));
    ALGEBRAIC(:,50) = ALGEBRAIC(:,46)./1000.00;
    ALGEBRAIC(:,52) = ALGEBRAIC(:,42);
    ALGEBRAIC(:,49) = ( CONSTANTS(:,48).*power(STATES(:,2), 2.00000))./(power(CONSTANTS(:,49), 2.00000)+power(STATES(:,2), 2.00000));
    ALGEBRAIC(:,51) =  (( (( (( CONSTANTS(:,50).*1.00000)./(power(CONSTANTS(:,51), 3.00000)+power(CONSTANTS(:,8), 3.00000))).*1.00000)./(CONSTANTS(:,52)+CONSTANTS(:,9))).*1.00000)./(1.00000+ CONSTANTS(:,53).*exp(( (CONSTANTS(:,54) - 1.00000).*ALGEBRAIC(:,72).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11))))).*( exp(( CONSTANTS(:,54).*ALGEBRAIC(:,72).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11))).*power(STATES(:,19), 3.00000).*CONSTANTS(:,9) -  exp(( (CONSTANTS(:,54) - 1.00000).*ALGEBRAIC(:,72).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11))).*power(CONSTANTS(:,8), 3.00000).*STATES(:,2));
    ALGEBRAIC(:,53) =  (( CONSTANTS(:,10).*CONSTANTS(:,11))./( 2.00000.*CONSTANTS(:,12))).*log(CONSTANTS(:,9)./STATES(:,2));
    ALGEBRAIC(:,55) =  CONSTANTS(:,55).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,53));
    ALGEBRAIC(:,12) = power(1.00000+( CONSTANTS(:,18).*CONSTANTS(:,20))./power(CONSTANTS(:,20)+STATES(:,2), 2.00000),  - 1.00000);
    ALGEBRAIC(:,45) = ( CONSTANTS(:,22).*STATES(:,2).*(CONSTANTS(:,33) - STATES(:,8))+ CONSTANTS(:,24).*STATES(:,2).*(CONSTANTS(:,32) - STATES(:,7))) - ( CONSTANTS(:,23).*STATES(:,8)+ CONSTANTS(:,25).*STATES(:,7));
    ALGEBRAIC(:,54) = ALGEBRAIC(:,46)./95000.0;
    ALGEBRAIC(:,56) = ALGEBRAIC(:,42)./50.0000;
    ALGEBRAIC(:,57) =  (( CONSTANTS(:,10).*CONSTANTS(:,11))./CONSTANTS(:,12)).*log(( 0.900000.*CONSTANTS(:,8)+ 0.100000.*CONSTANTS(:,7))./( 0.900000.*STATES(:,19)+ 0.100000.*STATES(:,28)));
    ALGEBRAIC(:,58) =  CONSTANTS(:,56).*STATES(:,20).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,57));
    ALGEBRAIC(:,59) =  CONSTANTS(:,57).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,57));
    ALGEBRAIC(:,68) = 1.00000./(1.00000+ 0.124500.*exp((  - 0.100000.*ALGEBRAIC(:,72).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11)))+ 0.0365000.*CONSTANTS(:,72).*exp((  - ALGEBRAIC(:,72).*CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11))));
    ALGEBRAIC(:,69) = ( (( CONSTANTS(:,66).*ALGEBRAIC(:,68).*1.00000)./(1.00000+power(CONSTANTS(:,67)./STATES(:,19), 1.50000))).*CONSTANTS(:,7))./(CONSTANTS(:,7)+CONSTANTS(:,68));
    ALGEBRAIC(:,60) =  (( CONSTANTS(:,10).*CONSTANTS(:,11))./CONSTANTS(:,12)).*log(CONSTANTS(:,7)./STATES(:,28));
    ALGEBRAIC(:,61) =  CONSTANTS(:,58).*power(STATES(:,29), 3.00000).*STATES(:,30).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60));
    ALGEBRAIC(:,62) =  CONSTANTS(:,59).*STATES(:,31).*STATES(:,32).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60));
    ALGEBRAIC(:,63) = ( (( 0.293800.*CONSTANTS(:,7))./(CONSTANTS(:,7)+210.000)).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60)))./(1.00000+exp( 0.0896000.*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60))));
    ALGEBRAIC(:,64) =  CONSTANTS(:,60).*power(STATES(:,33), 2.00000).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60));
    ALGEBRAIC(:,65) =  CONSTANTS(:,61).*STATES(:,34).*STATES(:,35).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60));
    ALGEBRAIC(:,66) =  CONSTANTS(:,62).*STATES(:,36).*STATES(:,37).*(ALGEBRAIC(:,72) - ALGEBRAIC(:,60));
    ALGEBRAIC(:,67) =  CONSTANTS(:,63).*STATES(:,38).*(ALGEBRAIC(:,72) -  (( CONSTANTS(:,10).*CONSTANTS(:,11))./CONSTANTS(:,12)).*log(( 0.980000.*CONSTANTS(:,7)+ 0.0200000.*CONSTANTS(:,8))./( 0.980000.*STATES(:,28)+ 0.0200000.*STATES(:,19))));
    ALGEBRAIC(:,1) = piecewise({t>=CONSTANTS(:,13)&t<=CONSTANTS(:,14)&(t - CONSTANTS(:,13)) -  floor((t - CONSTANTS(:,13))./CONSTANTS(:,15)).*CONSTANTS(:,15)<=CONSTANTS(:,16), CONSTANTS(:,17) }, 0.00000);
    ALGEBRAIC(:,70) = 0.200000./(1.00000+exp( - (ALGEBRAIC(:,72) - 46.7000)./7.80000));
    ALGEBRAIC(:,71) =  (( CONSTANTS(:,69).*ALGEBRAIC(:,70).*STATES(:,2))./(STATES(:,2)+CONSTANTS(:,71))).*(ALGEBRAIC(:,72) - CONSTANTS(:,70));
end

% Compute result of a piecewise function (for stimulation current)
function x = piecewise(cases, default)
    set = [0];
    for i = 1:2:length(cases)
        if (length(cases{i+1}) == 1)
            x(cases{i} & ~set,:) = cases{i+1};
        else
            x(cases{i} & ~set,:) = cases{i+1}(cases{i} & ~set);
        end
        set = set | cases{i};
        if(set), break, end
    end
    if (length(default) == 1)
        x(~set,:) = default;
    else
        x(~set,:) = default(~set);
    end
end

function VC = volt_clamp(t, holding_p, holding_t, P1, P1_t, P2)
    if t < holding_t
        VC = holding_p;
    elseif (t >= holding_t) && (t <= P1_t) 
        VC = P1;
    else
        VC = P2;
    end
end

function [algebraicVariableCount] = getAlgebraicVariableCount() 
    % Used later when setting a global variable with the number of algebraic variables.
    % There are a total of 41 entries in each of the rate and state variable arrays.
    % There are a total of 73 entries in the constant variable array.
    algebraicVariableCount =72;
end

function [STATES, CONSTANTS] = initConsts()
    CONSTANTS = []; STATES = [];

    % Table 8. Initial conditions
    STATES(:,1) = -82.4202;  % V; Membrane potential:mV
    STATES(:,2) = 0.115001;  % Cai; Myoplasmic Ca+ concentration:uM
    STATES(:,3) = 0.115001;  % Cass; Subspace SR Ca+ concentration:uM
    STATES(:,4) = 1299.5;  % CaJSR; JSR Ca+ concentration:uM
    STATES(:,5) = 1299.5;  % CaNSR; NSR Ca+ concentration:uM
    STATES(:,6) = 0;  % PRyR; RyR modulation factor
    STATES(:,7) = 11.2684;  % LTRPNCa; Concentration Ca+ bound low-affinity troponin-binding sites:uM
    STATES(:,8) = 125.29;  % HTEPNCa; Concentration Ca+ bound high-affinity troponin-binding sites:uM
    STATES(:,9) = 0.149102e-4;  % PO1; %Fraction of RyR channels in sate POl
    STATES(:,10) = 0.951726e-10;  % PO2; %Fraction of RyR channels in sate POl
    STATES(:,11) = 0.16774e-3;  % PC2; Fraction of RyR channels in sate PC2
    STATES(:,12) = 0.930308e-18;  % O; L-type Ca+ channel conducting state
    STATES(:,13) = 0.124216e-3;  % C2; L-type Ca+ channel closed state
    STATES(:,14) = 0.578679e-8;  % C3; L-type Ca+ channel closed state
    STATES(:,15) = 0.119816e-12;  % C4; L-type Ca+ channel closed state
    STATES(:,16) = 0.497923e-18;  % I1; L-type Ca+ channel inacticated state
    STATES(:,17) = 0.345847e-13;  % I2; L-type Ca+ channel inacticated state
    STATES(:,18) = 0.185106e-13;  % I3; L-type Ca+ channel inacticated state
    STATES(:,19) = 14237.1;  % Nai; Myoplasmic Na+ concentration:uM
    STATES(:,20) = 0.713483e-6;  % ONa; Open state of fast Na+ channel
    STATES(:,21) = 0.279132e-3;  % CNa1; Closed state of fast Na+ channel
    STATES(:,22) = 0.020752;  % CNa2; Closed state of fast Na+ channel
    STATES(:,23) = 0.673345e-6;  % I1Na; Slow inactivated state 1 of fast Na+ channel
    STATES(:,24) = 0.155787e-8;  % I2Na; Slow inactivated state 2 of fast Na+ channel
    STATES(:,25) = 0.153176e-3;  % IFNa; Fast inactivated state of fast Na+ channel
    STATES(:,26) = 0.0113879;  % ICNa2; Cloesd-inactivated state of fast Na+ channel
    STATES(:,27) = 0.34278;  % ICNa3; Cloesd-inactivated state of fast Na+ channel
    STATES(:,28) = 143720;  % Ki; Myoplasmic K+ concentration:uM
    STATES(:,29) = 0.265563e-2;  % ato_f; Gating variable for transient outward K+ current
    STATES(:,30) = 0.999977;  % ito_f; Gating variable for transient outward K+ current
    STATES(:,31) = 0.417069e-3;  % ato_s; Gating variable for transient outward K+ current
    STATES(:,32) = 0.998543;  % ito_s; Gating variable for transient outward K+ current
    STATES(:,33) = 0.262753e-3;  % nKs; Gating variable for slow delayed-rectifier K+ current
    STATES(:,34) = 0.417069e-3;  % aur; Gating variable for ultrarapidly activating delayed-rectifier K+ current
    STATES(:,35) = 0.998543;  % iur; Gating variable for ultrarapidly activating delayed-rectifier K+ current
    STATES(:,36) = 0.417069e-3;  % aKss; Gating variable for noninactivating steady-state K+ current
    STATES(:,37) = 1;  % iKss; Gating variable for noninactivating steady-state K+ current
    STATES(:,38) = 0.175298e-3;  % OK; mERG channel open state
    STATES(:,39) = 0.992513e-3;  % CK1; mERG channel closed state
    STATES(:,40) = 0.641229e-3;  % CK2; mERG channel closed state
    STATES(:,41) = 0.319129e-4;  % IK; mERG channel inactivated state
    
    % Table 2. Cell geometry parameters
    CONSTANTS(:,2) = 25.84e-6;  % Vmyo; Myoplasmic volume:ul
    CONSTANTS(:,3) = 0.12e-6;  % VJS; Juncional SR volume:ul
    CONSTANTS(:,4) = 2.098e-6;  % VNSR; Network SR volume:ul
    CONSTANTS(:,5) = 1.485e-9;  % Vss; Subspace volume:ul
    CONSTANTS(:,6) = 1.534e-4;  % Acap; Capacitive membrane area:cm^2
    
    % Table 3. Extracellular ion concentrations
    CONSTANTS(:,7) = 5400;  % Ko; Exracellular K+ concentration:uM
    CONSTANTS(:,8) = 140000;  % Nao; Exracellular Na+ concentration:uM
    CONSTANTS(:,9) = 1800;  % Cao; Exracellular Ca+ concentration:uM
    
    % Table 4. SR parameters
    CONSTANTS(:,26) = 4.5;  % nu1; Maximum RyR channel Ca+ permeability:ms^-1
    CONSTANTS(:,28) = 1.74e-5;  % nu2; Ca+ leak rate constant from the NSR:ms^-1
    CONSTANTS(:,27) = 20;  % tautr; Time constant for transfer from NSR to JSR:ms
    CONSTANTS(:,29) = 8;  % tauxfer; Time constant for transfer from subspace to myoplasm:ms
    CONSTANTS(:,30) = 0.45;  % nu3; SR Ca+-ATPase maximum pump rate:uM/ms
    CONSTANTS(:,31) = 0.5;  % Km_up; Half-saturation constant for SR Ca+-ATPase pump:uM
    CONSTANTS(:,35) = 0.006075;  % ka1; RyR Pc1-Po1 rate constant:uM^-4/ms
    CONSTANTS(:,36) = 0.07125;  % ka2; RyR Po1-Pc1 rate constant:ms^-1
    CONSTANTS(:,37) = 0.00405;  % kb1; RyR Po1-Po2 rate constant:uM^-3/ms
    CONSTANTS(:,38) = 0.965;  % kb2; RyR Po2-Po1 rate constant:ms^-1
    CONSTANTS(:,39) = 0.009;  % kc1; RyR Po1-Pc2 rate constant:ms^-1
    CONSTANTS(:,40) = 0.0008;  % kc2; RyR Pc2-Po1 rate constant:ms^-1
    CONSTANTS(:,41) = 3;  % m; RyR Ca+ cooperativity parameter Po1-Po2
    CONSTANTS(:,42) = 4;  % n; RyR Ca+ cooperativity parameter Pc1-Po1
    
    % Table 5. L-type Ca+ channel parameters
    CONSTANTS(:,44) = 0.1729;  % GCaL; Specific maximum conductivity for L-type Ca+ channel:mS/uF
    CONSTANTS(:,43) = 63;  % ECa; Reversal potential for L-type Ca+ channel:mV
    CONSTANTS(:,46) = 0.23324;  % Kpc_max; Maximum time constant for Ca+-induced inactivation:ms^-1
    CONSTANTS(:,47) = 20;  % Kpc_half; Half-saturation constant for Ca+-induced inactivation:uM
    CONSTANTS(:,45) = 0.0005;  % Kpcb; Voltage-insensitive rate constant for inactivation:ms^-1
    CONSTANTS(:,34) = 7;  % ICaL_max; Normalization constant for L-type Ca+ current:pA/pF 
    
    % Table 6. Buffering parameters
    CONSTANTS(:,32) = 70;  % LTRPN_tot; Total myoplasmic troponin low-affinity site concentration:uM
    CONSTANTS(:,33) = 140;  % HTRPN_tot; Total myoplasmic troponin high-affinity site concentration:uM
    CONSTANTS(:,22) = 0.00237;  % khtrpn1; Ca+ on rate constant for troponin high-affinity sites:uM^-1/ms
    CONSTANTS(:,23) = 3.2e-5;  % khtrpn2; Ca+ off rate constant for troponin high-affinity sites:ms^-1
    CONSTANTS(:,24) = 0.0327;  % kltrpn1; Ca+ on rate constant for troponin low-affinity sites:uM^-1/ms
    CONSTANTS(:,25) = 0.0196;  % kltrpn2; Ca+ off rate constant for troponin low-affinity sites:ms^-1
    CONSTANTS(:,18) = 50;  % CMDN_tot; Total myoplasmic calmodulin concentration:uM
    CONSTANTS(:,19) = 15000;  % CSQN_tot; Total junctional SR calsequestrin concentration:uM
    CONSTANTS(:,20) = 0.238;  % Km_CMDN; Ca+ half-saturaion constant for calmodulin:uM
    CONSTANTS(:,21) = 800;  % Km_CSQN; Ca+ half-saturaion constant for calsequestrin:uM
    
    % Table 7. Membrane current parameters
    CONSTANTS(:,12) = 96.5;  % F; Faraday constant:C/mmol
    CONSTANTS(:,11) = 298;  % T; Absolute temperature:K
    CONSTANTS(:,10) = 8.314;  % R; Ideal gas constant:J*mol^-1*K^-1
    CONSTANTS(:,50) = 292.8;  % kNaCa; Scaling factor of Na+/Ca+ exchange:pA/pF
    CONSTANTS(:,51) = 87500;  % Km_Na; Na+ half-saturation constant for Na+/Ca+ exchange:uM
    CONSTANTS(:,52) = 1380;  %  Km_Ca; Ca+ half-saturation constant for Na+/Ca+ exchange:uM
    CONSTANTS(:,53) = 0.1;  % ksat; Na+/Ca+ exchange saaturation factor at very negative potentials
    CONSTANTS(:,54) = 0.35;  % yita; Contals voltage dependence of Na+/Ca+ exchange
    CONSTANTS(:,66) = 0.88;  % INaK_max; Maximum Na+/K+ exchange current:pA/pF 
    CONSTANTS(:,67) = 21000;  % Km_Nai; Na+ half-saturation constant for Na+/K+ exchange currant:uM
    CONSTANTS(:,68) = 1500;  % Km_Ko; K+ half-saturation constant for Na+/K+ exchange currant:uM
    CONSTANTS(:,48) = 1;  % Ip_Ca_max; Maximun Ca+ pump current:pA/pF
    CONSTANTS(:,49) = 0.5;  % Km_p_Ca; Ca+ half-saturation constant for Ca+pump current:uM
    CONSTANTS(:,55) = 0.000367;  % GCab; Maximun background Ca+ current conductance:mS/uF
    CONSTANTS(:,56) = 13;  % GNa; Maximun fast Na+ current conductance:mS/uF
    CONSTANTS(:,57) = 0.0026;  % GNab; Maximun background Na+ current conductance:mS/uF
    CONSTANTS(:,60) = 0.00575;  % GKs; Maximum slow delayed-rectifier K+ current conductance:mS/uF
    
    % Table 7. Parameters for apex
    CONSTANTS(:,58) = 0.4067;  % GKtof; Maximum transient outward K+ current conductance(apex):mS/uF
    CONSTANTS(:,59) = 0;  % GKtos; Maximum transient outward K+ current conductance(apex):mS/uF
    CONSTANTS(:,61) = 0.16;  % GKur; Maximum ultrarapidly delayed-rectifier K+ current conductance(apex):mS/uF
    CONSTANTS(:,62) = 0.05;  % GKss; Maximum noninactivating steady-state K+ current conductance(apex):mS/uF
    
    % Table 7. Membrane current parameters
    CONSTANTS(:,63) = 0.078;  % GKr; Maximum rapid delayed-rectifier K+ current conductance:mS/uF
    CONSTANTS(:,65) = 0.023761;  % kf; Rate constant for rapid delayed-rectifier K+ current:ms^-1
    CONSTANTS(:,64) = 0.036778;  % kb; Rate constant for rapid delayed-rectifier K+ current:ms^-1
    CONSTANTS(:,69) = 10;  % GCl_Ca; Maximum Ca+-activated Cl- current conductance:mS/uF
    CONSTANTS(:,71) = 10;  % Km_Cl; Half-saturaon constant for Ca+-activated Cl- current:uM
    CONSTANTS(:,70) = -40;  % ECl; Reversal potential for Ca+-activated Cl- current:mV
    
    % Parameters not in the Simulink model
    CONSTANTS(:,1) = 1;
    CONSTANTS(:,13) = 20;
    CONSTANTS(:,14) = 100000;
    CONSTANTS(:,15) = 71.43;
    CONSTANTS(:,16) = 0.5;
    CONSTANTS(:,17) = -80;
    CONSTANTS(:,72) =  (1.00000./7.00000).*(exp(CONSTANTS(:,8)./67300.0) - 1.00000);
    CONSTANTS(:,72) = 0.00000;

    if (isempty(STATES)), warning('Initial values for states not set'); end
end
