
function [VOI, STATES, ALGEBRAIC, CONSTANTS] = mainFunction()
    % This is the "main function".  In Matlab, things work best if you rename this function to match the filename.
   [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel();
end

function [algebraicVariableCount] = getAlgebraicVariableCount() 
    % Used later when setting a global variable with the number of algebraic variables.
    % Note: This is not the "main method".  
    algebraicVariableCount =54;
end
% There are a total of 26 entries in each of the rate and state variable arrays.
% There are a total of 70 entries in the constant variable array.
%

function [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel()
    % Create ALGEBRAIC of correct size
    global algebraicVariableCount;  algebraicVariableCount = getAlgebraicVariableCount();
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts;

    % Set timespan to solve over 
    tspan = [0, 10];

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 1);

    % Solve model with ODE solver
    [VOI, STATES] = ode15s(@(VOI, STATES)computeRates(VOI, STATES, CONSTANTS), tspan, INIT_STATES, options);

    % Compute algebraic variables
    [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI);

    % Plot state variables against variable of integration
    [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends();
    figure();
    plot(VOI, STATES);
    xlabel(LEGEND_VOI);
    l = legend(LEGEND_STATES);
    set(l,'Interpreter','none');
end

function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends()
    LEGEND_STATES = ''; LEGEND_ALGEBRAIC = ''; LEGEND_VOI = ''; LEGEND_CONSTANTS = '';
    LEGEND_VOI = strpad('time in component environment (second)');
    LEGEND_STATES(:,1) = strpad('V in component membrane (millivolt)');
    LEGEND_CONSTANTS(:,1) = strpad('R in component membrane (millijoule_per_mole_kelvin)');
    LEGEND_CONSTANTS(:,2) = strpad('T in component membrane (kelvin)');
    LEGEND_CONSTANTS(:,3) = strpad('F in component membrane (coulomb_per_mole)');
    LEGEND_CONSTANTS(:,4) = strpad('Cm in component membrane (microF)');
    LEGEND_ALGEBRAIC(:,26) = strpad('i_Na in component sodium_current (nanoA)');
    LEGEND_ALGEBRAIC(:,27) = strpad('i_Ca_L in component L_type_Ca_channel (nanoA)');
    LEGEND_ALGEBRAIC(:,29) = strpad('i_t in component Ca_independent_transient_outward_K_current (nanoA)');
    LEGEND_ALGEBRAIC(:,30) = strpad('i_ss in component steady_state_outward_K_current (nanoA)');
    LEGEND_ALGEBRAIC(:,35) = strpad('i_f in component hyperpolarisation_activated_current (nanoA)');
    LEGEND_ALGEBRAIC(:,32) = strpad('i_K1 in component inward_rectifier (nanoA)');
    LEGEND_ALGEBRAIC(:,39) = strpad('i_B in component background_currents (nanoA)');
    LEGEND_ALGEBRAIC(:,40) = strpad('i_NaK in component sodium_potassium_pump (nanoA)');
    LEGEND_ALGEBRAIC(:,42) = strpad('i_NaCa in component Na_Ca_ion_exchanger_current (nanoA)');
    LEGEND_ALGEBRAIC(:,41) = strpad('i_Ca_P in component sarcolemmal_calcium_pump_current (nanoA)');
    LEGEND_ALGEBRAIC(:,14) = strpad('i_Stim in component membrane (nanoA)');
    LEGEND_ALGEBRAIC(:,31) = strpad('i_K_slow in component slowly_inactivating_delayed_rectifier_K_current (nanoA)');
    LEGEND_CONSTANTS(:,5) = strpad('stim_period in component membrane (second)');
    LEGEND_CONSTANTS(:,6) = strpad('stim_duration in component membrane (second)');
    LEGEND_CONSTANTS(:,7) = strpad('stim_amplitude in component membrane (nanoA)');
    LEGEND_ALGEBRAIC(:,25) = strpad('E_Na in component sodium_current (millivolt)');
    LEGEND_CONSTANTS(:,8) = strpad('g_Na in component sodium_current (microS)');
    LEGEND_CONSTANTS(:,66) = strpad('g_Na_endo in component sodium_current (microS)');
    LEGEND_STATES(:,2) = strpad('Na_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_CONSTANTS(:,9) = strpad('Na_o in component standard_ionic_concentrations (millimolar)');
    LEGEND_STATES(:,3) = strpad('m in component sodium_current_m_gate (dimensionless)');
    LEGEND_STATES(:,4) = strpad('h in component sodium_current_h_gate (dimensionless)');
    LEGEND_STATES(:,5) = strpad('j in component sodium_current_j_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,1) = strpad('m_infinity in component sodium_current_m_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,15) = strpad('tau_m in component sodium_current_m_gate (second)');
    LEGEND_ALGEBRAIC(:,2) = strpad('h_infinity in component sodium_current_h_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,16) = strpad('tau_h in component sodium_current_h_gate (second)');
    LEGEND_ALGEBRAIC(:,3) = strpad('j_infinity in component sodium_current_j_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,17) = strpad('tau_j in component sodium_current_j_gate (second)');
    LEGEND_CONSTANTS(:,10) = strpad('g_Ca_L in component L_type_Ca_channel (microS)');
    LEGEND_CONSTANTS(:,11) = strpad('E_Ca_L in component L_type_Ca_channel (millivolt)');
    LEGEND_STATES(:,6) = strpad('Ca_ss in component intracellular_ion_concentrations (millimolar)');
    LEGEND_STATES(:,7) = strpad('d in component L_type_Ca_channel_d_gate (dimensionless)');
    LEGEND_STATES(:,8) = strpad('f_11 in component L_type_Ca_channel_f_11_gate (dimensionless)');
    LEGEND_STATES(:,9) = strpad('f_12 in component L_type_Ca_channel_f_12_gate (dimensionless)');
    LEGEND_STATES(:,10) = strpad('Ca_inact in component L_type_Ca_channel_Ca_inact_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,4) = strpad('d_infinity in component L_type_Ca_channel_d_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,18) = strpad('tau_d in component L_type_Ca_channel_d_gate (second)');
    LEGEND_ALGEBRAIC(:,5) = strpad('f_11_infinity in component L_type_Ca_channel_f_11_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,19) = strpad('tau_f_11 in component L_type_Ca_channel_f_11_gate (second)');
    LEGEND_ALGEBRAIC(:,6) = strpad('f_12_infinity in component L_type_Ca_channel_f_12_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,20) = strpad('tau_f_12 in component L_type_Ca_channel_f_12_gate (second)');
    LEGEND_CONSTANTS(:,12) = strpad('tau_Ca_inact in component L_type_Ca_channel_Ca_inact_gate (second)');
    LEGEND_ALGEBRAIC(:,7) = strpad('Ca_inact_infinity in component L_type_Ca_channel_Ca_inact_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,28) = strpad('E_K in component Ca_independent_transient_outward_K_current (millivolt)');
    LEGEND_CONSTANTS(:,13) = strpad('g_t in component Ca_independent_transient_outward_K_current (microS)');
    LEGEND_CONSTANTS(:,14) = strpad('K_o in component standard_ionic_concentrations (millimolar)');
    LEGEND_STATES(:,11) = strpad('K_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_STATES(:,12) = strpad('r in component Ca_independent_transient_outward_K_current_r_gate (dimensionless)');
    LEGEND_STATES(:,13) = strpad('s in component Ca_independent_transient_outward_K_current_s_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,21) = strpad('tau_r in component Ca_independent_transient_outward_K_current_r_gate (second)');
    LEGEND_ALGEBRAIC(:,8) = strpad('r_infinity in component Ca_independent_transient_outward_K_current_r_gate (dimensionless)');
    LEGEND_CONSTANTS(:,67) = strpad('tau_s in component Ca_independent_transient_outward_K_current_s_gate (second)');
    LEGEND_ALGEBRAIC(:,9) = strpad('s_infinity in component Ca_independent_transient_outward_K_current_s_gate (dimensionless)');
    LEGEND_CONSTANTS(:,15) = strpad('g_ss in component steady_state_outward_K_current (microS)');
    LEGEND_STATES(:,14) = strpad('r_ss in component steady_state_outward_K_current_r_ss_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,22) = strpad('tau_r_ss in component steady_state_outward_K_current_r_ss_gate (second)');
    LEGEND_ALGEBRAIC(:,10) = strpad('r_ss_infinity in component steady_state_outward_K_current_r_ss_gate (dimensionless)');
    LEGEND_CONSTANTS(:,16) = strpad('g_K_slow in component slowly_inactivating_delayed_rectifier_K_current (microS)');
    LEGEND_STATES(:,15) = strpad('r_K_slow in component slowly_inactivating_delayed_rectifier_K_current_r_K_slow_gate (dimensionless)');
    LEGEND_STATES(:,16) = strpad('s_K_slow in component slowly_inactivating_delayed_rectifier_K_current_s_K_slow_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,23) = strpad('tau_r_K_slow in component slowly_inactivating_delayed_rectifier_K_current_r_K_slow_gate (second)');
    LEGEND_ALGEBRAIC(:,11) = strpad('r_K_slow_infinity in component slowly_inactivating_delayed_rectifier_K_current_r_K_slow_gate (dimensionless)');
    LEGEND_CONSTANTS(:,68) = strpad('tau_s_K_slow in component slowly_inactivating_delayed_rectifier_K_current_s_K_slow_gate (second)');
    LEGEND_ALGEBRAIC(:,12) = strpad('s_K_slow_infinity in component slowly_inactivating_delayed_rectifier_K_current_s_K_slow_gate (dimensionless)');
    LEGEND_CONSTANTS(:,17) = strpad('g_K1 in component inward_rectifier (microS)');
    LEGEND_ALGEBRAIC(:,33) = strpad('i_f_Na in component hyperpolarisation_activated_current (nanoA)');
    LEGEND_ALGEBRAIC(:,34) = strpad('i_f_K in component hyperpolarisation_activated_current (nanoA)');
    LEGEND_CONSTANTS(:,18) = strpad('g_f in component hyperpolarisation_activated_current (microS)');
    LEGEND_CONSTANTS(:,19) = strpad('f_Na in component hyperpolarisation_activated_current (dimensionless)');
    LEGEND_CONSTANTS(:,69) = strpad('f_K in component hyperpolarisation_activated_current (dimensionless)');
    LEGEND_STATES(:,17) = strpad('y in component hyperpolarisation_activated_current_y_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,24) = strpad('tau_y in component hyperpolarisation_activated_current_y_gate (second)');
    LEGEND_ALGEBRAIC(:,13) = strpad('y_infinity in component hyperpolarisation_activated_current_y_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,36) = strpad('i_B_Na in component background_currents (nanoA)');
    LEGEND_ALGEBRAIC(:,37) = strpad('i_B_Ca in component background_currents (nanoA)');
    LEGEND_ALGEBRAIC(:,38) = strpad('i_B_K in component background_currents (nanoA)');
    LEGEND_CONSTANTS(:,20) = strpad('g_B_Na in component background_currents (microS)');
    LEGEND_CONSTANTS(:,21) = strpad('g_B_Ca in component background_currents (microS)');
    LEGEND_CONSTANTS(:,22) = strpad('g_B_K in component background_currents (microS)');
    LEGEND_CONSTANTS(:,23) = strpad('E_Ca in component background_currents (millivolt)');
    LEGEND_CONSTANTS(:,24) = strpad('Ca_o in component standard_ionic_concentrations (millimolar)');
    LEGEND_STATES(:,18) = strpad('Ca_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_CONSTANTS(:,25) = strpad('i_NaK_max in component sodium_potassium_pump (nanoA)');
    LEGEND_CONSTANTS(:,26) = strpad('K_m_K in component sodium_potassium_pump (millimolar)');
    LEGEND_CONSTANTS(:,27) = strpad('K_m_Na in component sodium_potassium_pump (millimolar)');
    LEGEND_CONSTANTS(:,70) = strpad('sigma in component sodium_potassium_pump (dimensionless)');
    LEGEND_CONSTANTS(:,28) = strpad('i_Ca_P_max in component sarcolemmal_calcium_pump_current (nanoA)');
    LEGEND_CONSTANTS(:,29) = strpad('K_NaCa in component Na_Ca_ion_exchanger_current (millimolar_4)');
    LEGEND_CONSTANTS(:,30) = strpad('d_NaCa in component Na_Ca_ion_exchanger_current (millimolar_4)');
    LEGEND_CONSTANTS(:,31) = strpad('gamma_NaCa in component Na_Ca_ion_exchanger_current (dimensionless)');
    LEGEND_ALGEBRAIC(:,43) = strpad('J_rel in component SR_Ca_release_channel (millimolar_per_second)');
    LEGEND_CONSTANTS(:,32) = strpad('v1 in component SR_Ca_release_channel (per_second)');
    LEGEND_CONSTANTS(:,33) = strpad('k_a_plus in component SR_Ca_release_channel (millimolar4_per_second)');
    LEGEND_CONSTANTS(:,34) = strpad('k_a_minus in component SR_Ca_release_channel (per_second)');
    LEGEND_CONSTANTS(:,35) = strpad('k_b_plus in component SR_Ca_release_channel (millimolar3_per_second)');
    LEGEND_CONSTANTS(:,36) = strpad('k_b_minus in component SR_Ca_release_channel (per_second)');
    LEGEND_CONSTANTS(:,37) = strpad('k_c_plus in component SR_Ca_release_channel (per_second)');
    LEGEND_CONSTANTS(:,38) = strpad('k_c_minus in component SR_Ca_release_channel (per_second)');
    LEGEND_STATES(:,19) = strpad('P_O1 in component SR_Ca_release_channel (dimensionless)');
    LEGEND_STATES(:,20) = strpad('P_O2 in component SR_Ca_release_channel (dimensionless)');
    LEGEND_STATES(:,21) = strpad('P_C1 in component SR_Ca_release_channel (dimensionless)');
    LEGEND_STATES(:,22) = strpad('P_C2 in component SR_Ca_release_channel (dimensionless)');
    LEGEND_CONSTANTS(:,39) = strpad('n in component SR_Ca_release_channel (dimensionless)');
    LEGEND_CONSTANTS(:,40) = strpad('m in component SR_Ca_release_channel (dimensionless)');
    LEGEND_STATES(:,23) = strpad('Ca_JSR in component intracellular_ion_concentrations (millimolar)');
    LEGEND_ALGEBRAIC(:,46) = strpad('J_up in component SERCA2a_pump (millimolar_per_second)');
    LEGEND_CONSTANTS(:,41) = strpad('K_fb in component SERCA2a_pump (millimolar)');
    LEGEND_CONSTANTS(:,42) = strpad('K_rb in component SERCA2a_pump (millimolar)');
    LEGEND_ALGEBRAIC(:,44) = strpad('fb in component SERCA2a_pump (dimensionless)');
    LEGEND_ALGEBRAIC(:,45) = strpad('rb in component SERCA2a_pump (dimensionless)');
    LEGEND_CONSTANTS(:,43) = strpad('Vmaxf in component SERCA2a_pump (millimolar_per_second)');
    LEGEND_CONSTANTS(:,44) = strpad('Vmaxr in component SERCA2a_pump (millimolar_per_second)');
    LEGEND_CONSTANTS(:,45) = strpad('K_SR in component SERCA2a_pump (dimensionless)');
    LEGEND_CONSTANTS(:,46) = strpad('N_fb in component SERCA2a_pump (dimensionless)');
    LEGEND_CONSTANTS(:,47) = strpad('N_rb in component SERCA2a_pump (dimensionless)');
    LEGEND_STATES(:,24) = strpad('Ca_NSR in component intracellular_ion_concentrations (millimolar)');
    LEGEND_ALGEBRAIC(:,48) = strpad('J_tr in component intracellular_and_SR_Ca_fluxes (millimolar_per_second)');
    LEGEND_ALGEBRAIC(:,47) = strpad('J_xfer in component intracellular_and_SR_Ca_fluxes (millimolar_per_second)');
    LEGEND_ALGEBRAIC(:,53) = strpad('J_trpn in component intracellular_and_SR_Ca_fluxes (millimolar_per_second)');
    LEGEND_CONSTANTS(:,48) = strpad('tau_tr in component intracellular_and_SR_Ca_fluxes (second)');
    LEGEND_CONSTANTS(:,49) = strpad('tau_xfer in component intracellular_and_SR_Ca_fluxes (second)');
    LEGEND_STATES(:,25) = strpad('HTRPNCa in component intracellular_and_SR_Ca_fluxes (millimolar)');
    LEGEND_STATES(:,26) = strpad('LTRPNCa in component intracellular_and_SR_Ca_fluxes (millimolar)');
    LEGEND_ALGEBRAIC(:,49) = strpad('J_HTRPNCa in component intracellular_and_SR_Ca_fluxes (millimolar_per_second)');
    LEGEND_ALGEBRAIC(:,52) = strpad('J_LTRPNCa in component intracellular_and_SR_Ca_fluxes (millimolar_per_second)');
    LEGEND_CONSTANTS(:,50) = strpad('HTRPN_tot in component intracellular_and_SR_Ca_fluxes (millimolar)');
    LEGEND_CONSTANTS(:,51) = strpad('LTRPN_tot in component intracellular_and_SR_Ca_fluxes (millimolar)');
    LEGEND_CONSTANTS(:,52) = strpad('k_htrpn_plus in component intracellular_and_SR_Ca_fluxes (millimolar_per_second)');
    LEGEND_CONSTANTS(:,53) = strpad('k_htrpn_minus in component intracellular_and_SR_Ca_fluxes (per_second)');
    LEGEND_CONSTANTS(:,54) = strpad('k_ltrpn_plus in component intracellular_and_SR_Ca_fluxes (millimolar_per_second)');
    LEGEND_CONSTANTS(:,55) = strpad('k_ltrpn_minus in component intracellular_and_SR_Ca_fluxes (per_second)');
    LEGEND_CONSTANTS(:,56) = strpad('V_myo in component intracellular_ion_concentrations (micro_litre)');
    LEGEND_CONSTANTS(:,57) = strpad('V_JSR in component intracellular_ion_concentrations (micro_litre)');
    LEGEND_CONSTANTS(:,58) = strpad('V_NSR in component intracellular_ion_concentrations (micro_litre)');
    LEGEND_CONSTANTS(:,59) = strpad('V_SS in component intracellular_ion_concentrations (micro_litre)');
    LEGEND_CONSTANTS(:,60) = strpad('K_mCMDN in component intracellular_ion_concentrations (millimolar)');
    LEGEND_CONSTANTS(:,61) = strpad('K_mCSQN in component intracellular_ion_concentrations (millimolar)');
    LEGEND_CONSTANTS(:,62) = strpad('K_mEGTA in component intracellular_ion_concentrations (millimolar)');
    LEGEND_CONSTANTS(:,63) = strpad('CMDN_tot in component intracellular_ion_concentrations (millimolar)');
    LEGEND_CONSTANTS(:,64) = strpad('CSQN_tot in component intracellular_ion_concentrations (millimolar)');
    LEGEND_CONSTANTS(:,65) = strpad('EGTA_tot in component intracellular_ion_concentrations (millimolar)');
    LEGEND_ALGEBRAIC(:,54) = strpad('beta_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_ALGEBRAIC(:,50) = strpad('beta_SS in component intracellular_ion_concentrations (millimolar)');
    LEGEND_ALGEBRAIC(:,51) = strpad('beta_JSR in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES(:,1) = strpad('d/dt V in component membrane (millivolt)');
    LEGEND_RATES(:,3) = strpad('d/dt m in component sodium_current_m_gate (dimensionless)');
    LEGEND_RATES(:,4) = strpad('d/dt h in component sodium_current_h_gate (dimensionless)');
    LEGEND_RATES(:,5) = strpad('d/dt j in component sodium_current_j_gate (dimensionless)');
    LEGEND_RATES(:,7) = strpad('d/dt d in component L_type_Ca_channel_d_gate (dimensionless)');
    LEGEND_RATES(:,8) = strpad('d/dt f_11 in component L_type_Ca_channel_f_11_gate (dimensionless)');
    LEGEND_RATES(:,9) = strpad('d/dt f_12 in component L_type_Ca_channel_f_12_gate (dimensionless)');
    LEGEND_RATES(:,10) = strpad('d/dt Ca_inact in component L_type_Ca_channel_Ca_inact_gate (dimensionless)');
    LEGEND_RATES(:,12) = strpad('d/dt r in component Ca_independent_transient_outward_K_current_r_gate (dimensionless)');
    LEGEND_RATES(:,13) = strpad('d/dt s in component Ca_independent_transient_outward_K_current_s_gate (dimensionless)');
    LEGEND_RATES(:,14) = strpad('d/dt r_ss in component steady_state_outward_K_current_r_ss_gate (dimensionless)');
    LEGEND_RATES(:,15) = strpad('d/dt r_K_slow in component slowly_inactivating_delayed_rectifier_K_current_r_K_slow_gate (dimensionless)');
    LEGEND_RATES(:,16) = strpad('d/dt s_K_slow in component slowly_inactivating_delayed_rectifier_K_current_s_K_slow_gate (dimensionless)');
    LEGEND_RATES(:,17) = strpad('d/dt y in component hyperpolarisation_activated_current_y_gate (dimensionless)');
    LEGEND_RATES(:,21) = strpad('d/dt P_C1 in component SR_Ca_release_channel (dimensionless)');
    LEGEND_RATES(:,19) = strpad('d/dt P_O1 in component SR_Ca_release_channel (dimensionless)');
    LEGEND_RATES(:,20) = strpad('d/dt P_O2 in component SR_Ca_release_channel (dimensionless)');
    LEGEND_RATES(:,22) = strpad('d/dt P_C2 in component SR_Ca_release_channel (dimensionless)');
    LEGEND_RATES(:,25) = strpad('d/dt HTRPNCa in component intracellular_and_SR_Ca_fluxes (millimolar)');
    LEGEND_RATES(:,26) = strpad('d/dt LTRPNCa in component intracellular_and_SR_Ca_fluxes (millimolar)');
    LEGEND_RATES(:,18) = strpad('d/dt Ca_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES(:,2) = strpad('d/dt Na_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES(:,11) = strpad('d/dt K_i in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES(:,6) = strpad('d/dt Ca_ss in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES(:,23) = strpad('d/dt Ca_JSR in component intracellular_ion_concentrations (millimolar)');
    LEGEND_RATES(:,24) = strpad('d/dt Ca_NSR in component intracellular_ion_concentrations (millimolar)');
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end

function [STATES, CONSTANTS] = initConsts()
    VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];
    STATES(:,1) = -80.50146;
    CONSTANTS(:,1) = 8314.5;
    CONSTANTS(:,2) = 295;
    CONSTANTS(:,3) = 96487;
    CONSTANTS(:,4) = 0.0001;
    CONSTANTS(:,5) = 1;
    CONSTANTS(:,6) = 5e-3;
    CONSTANTS(:,7) = -0.6;
    CONSTANTS(:,8) = 1.064;
    STATES(:,2) = 10.73519;
    CONSTANTS(:,9) = 140;
    STATES(:,3) = 0.004164108;
    STATES(:,4) = 0.6735613;
    STATES(:,5) = 0.6729362;
    CONSTANTS(:,10) = 0.0341;
    CONSTANTS(:,11) = 65;
    STATES(:,6) = 0.00008737212;
    STATES(:,7) = 0.000002171081;
    STATES(:,8) = 0.9999529;
    STATES(:,9) = 0.9999529;
    STATES(:,10) = 0.9913102;
    CONSTANTS(:,12) = 0.009;
    CONSTANTS(:,13) = 0.033;
    CONSTANTS(:,14) = 5.4;
    STATES(:,11) = 139.2751;
    STATES(:,12) = 0.002191519;
    STATES(:,13) = 0.9842542;
    CONSTANTS(:,15) = 0.005;
    STATES(:,14) = 0.002907171;
    CONSTANTS(:,16) = 0.014;
    STATES(:,15) = 0.642;
    STATES(:,16) = 0.314;
    CONSTANTS(:,17) = 0.024;
    CONSTANTS(:,18) = 0.00145;
    CONSTANTS(:,19) = 0.2;
    STATES(:,17) = 0.003578708;
    CONSTANTS(:,20) = 0.00008015;
    CONSTANTS(:,21) = 0.0000324;
    CONSTANTS(:,22) = 0.000138;
    CONSTANTS(:,23) = 65;
    CONSTANTS(:,24) = 1.2;
    STATES(:,18) = 0.00007901351;
    CONSTANTS(:,25) = 0.08;
    CONSTANTS(:,26) = 1.5;
    CONSTANTS(:,27) = 10;
    CONSTANTS(:,28) = 0.004;
    CONSTANTS(:,29) = 0.000009984;
    CONSTANTS(:,30) = 0.0001;
    CONSTANTS(:,31) = 0.5;
    CONSTANTS(:,32) = 1.8e3;
    CONSTANTS(:,33) = 12.15e12;
    CONSTANTS(:,34) = 576;
    CONSTANTS(:,35) = 4.05e9;
    CONSTANTS(:,36) = 1930;
    CONSTANTS(:,37) = 100;
    CONSTANTS(:,38) = 0.8;
    STATES(:,19) = 0.0004327548;
    STATES(:,20) = 0.000000000606254;
    STATES(:,21) = 0.6348229;
    STATES(:,22) = 0.3647471;
    CONSTANTS(:,39) = 4;
    CONSTANTS(:,40) = 3;
    STATES(:,23) = 0.06607948;
    CONSTANTS(:,41) = 0.000168;
    CONSTANTS(:,42) = 3.29;
    CONSTANTS(:,43) = 0.04;
    CONSTANTS(:,44) = 0.9;
    CONSTANTS(:,45) = 1;
    CONSTANTS(:,46) = 1.2;
    CONSTANTS(:,47) = 1;
    STATES(:,24) = 0.06600742;
    CONSTANTS(:,48) = 0.0005747;
    CONSTANTS(:,49) = 0.0267;
    STATES(:,25) = 1.394301e-1;
    STATES(:,26) = 5.1619e-3;
    CONSTANTS(:,50) = 0.14;
    CONSTANTS(:,51) = 0.07;
    CONSTANTS(:,52) = 200000;
    CONSTANTS(:,53) = 0.066;
    CONSTANTS(:,54) = 40000;
    CONSTANTS(:,55) = 40;
    CONSTANTS(:,56) = 0.00000936;
    CONSTANTS(:,57) = 0.000000056;
    CONSTANTS(:,58) = 0.000000504;
    CONSTANTS(:,59) = 0.0000000012;
    CONSTANTS(:,60) = 0.00238;
    CONSTANTS(:,61) = 0.8;
    CONSTANTS(:,62) = 0.00015;
    CONSTANTS(:,63) = 0.05;
    CONSTANTS(:,64) = 15;
    CONSTANTS(:,65) = 0;
    CONSTANTS(:,66) =  1.33000.*CONSTANTS(:,8);
    CONSTANTS(:,67) = 0.0572000;
    CONSTANTS(:,68) = 1.17400;
    CONSTANTS(:,69) = 1.00000 - CONSTANTS(:,19);
    CONSTANTS(:,70) = (exp(CONSTANTS(:,9)./67.3000) - 1.00000)./7.00000;
    if (isempty(STATES)), warning('Initial values for states not set');, end
end

function [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS)
    global algebraicVariableCount;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        ALGEBRAIC = zeros(1, algebraicVariableCount);
        utilOnes = 1;
    else
        statesRowCount = statesSize(1);
        ALGEBRAIC = zeros(statesRowCount, algebraicVariableCount);
        RATES = zeros(statesRowCount, statesColumnCount);
        utilOnes = ones(statesRowCount, 1);
    end
    RATES(:,21) =   - CONSTANTS(:,33).*power(STATES(:,6), CONSTANTS(:,39)).*STATES(:,21)+ CONSTANTS(:,34).*STATES(:,19);
    RATES(:,19) = ( CONSTANTS(:,33).*power(STATES(:,6), CONSTANTS(:,39)).*STATES(:,21) - ( CONSTANTS(:,34).*STATES(:,19)+ CONSTANTS(:,35).*power(STATES(:,6), CONSTANTS(:,40)).*STATES(:,19)+ CONSTANTS(:,37).*STATES(:,19)))+ CONSTANTS(:,36).*STATES(:,20)+ CONSTANTS(:,38).*STATES(:,22);
    RATES(:,20) =  CONSTANTS(:,35).*power(STATES(:,6), CONSTANTS(:,40)).*STATES(:,19) -  CONSTANTS(:,36).*STATES(:,20);
    RATES(:,22) =  CONSTANTS(:,37).*STATES(:,19) -  CONSTANTS(:,38).*STATES(:,22);
    ALGEBRAIC(:,7) = 1.00000./(1.00000+STATES(:,6)./0.0100000);
    RATES(:,10) = (ALGEBRAIC(:,7) - STATES(:,10))./CONSTANTS(:,12);
    ALGEBRAIC(:,9) = 1.00000./(1.00000+exp((STATES(:,1)+24.8000)./3.50000));
    RATES(:,13) = (ALGEBRAIC(:,9) - STATES(:,13))./CONSTANTS(:,67);
    ALGEBRAIC(:,12) = 1.00000./(1.00000+exp((STATES(:,1)+37.6000)./5.90000));
    RATES(:,16) = (ALGEBRAIC(:,12) - STATES(:,16))./CONSTANTS(:,68);
    ALGEBRAIC(:,1) = 1.00000./(1.00000+exp((STATES(:,1)+45.0000)./ - 6.50000));
    ALGEBRAIC(:,15) = 0.00136000./(( 0.320000.*(STATES(:,1)+47.1300))./(1.00000 - exp(  - 0.100000.*(STATES(:,1)+47.1300)))+ 0.0800000.*exp( - STATES(:,1)./11.0000));
    RATES(:,3) = (ALGEBRAIC(:,1) - STATES(:,3))./ALGEBRAIC(:,15);
    ALGEBRAIC(:,2) = 1.00000./(1.00000+exp((STATES(:,1)+76.1000)./6.07000));
    ALGEBRAIC(:,16) = piecewise({STATES(:,1)>= - 40.0000,  0.000453700.*(1.00000+exp( - (STATES(:,1)+10.6600)./11.1000)) }, 0.00349000./( 0.135000.*exp( - (STATES(:,1)+80.0000)./6.80000)+ 3.56000.*exp( 0.0790000.*STATES(:,1))+ 310000..*exp( 0.350000.*STATES(:,1))));
    RATES(:,4) = (ALGEBRAIC(:,2) - STATES(:,4))./ALGEBRAIC(:,16);
    ALGEBRAIC(:,3) = 1.00000./(1.00000+exp((STATES(:,1)+76.1000)./6.07000));
    ALGEBRAIC(:,17) = piecewise({STATES(:,1)>= - 40.0000, ( 0.0116300.*(1.00000+exp(  - 0.100000.*(STATES(:,1)+32.0000))))./exp(  - 2.53500e-07.*STATES(:,1)) }, 0.00349000./( ((STATES(:,1)+37.7800)./(1.00000+exp( 0.311000.*(STATES(:,1)+79.2300)))).*(  - 127140..*exp( 0.244400.*STATES(:,1)) -  3.47400e-05.*exp(  - 0.0439100.*STATES(:,1)))+( 0.121200.*exp(  - 0.0105200.*STATES(:,1)))./(1.00000+exp(  - 0.137800.*(STATES(:,1)+40.1400)))));
    RATES(:,5) = (ALGEBRAIC(:,3) - STATES(:,5))./ALGEBRAIC(:,17);
    ALGEBRAIC(:,4) = 1.00000./(1.00000+exp((STATES(:,1)+15.3000)./ - 5.00000));
    ALGEBRAIC(:,18) =  0.00305000.*exp(  - 0.00450000.*power(STATES(:,1)+7.00000, 2.00000))+ 0.00105000.*exp(  - 0.00200000.*power(STATES(:,1) - 18.0000, 2.00000))+0.000250000;
    RATES(:,7) = (ALGEBRAIC(:,4) - STATES(:,7))./ALGEBRAIC(:,18);
    ALGEBRAIC(:,5) = 1.00000./(1.00000+exp((STATES(:,1)+26.7000)./5.40000));
    ALGEBRAIC(:,19) =  0.105000.*exp( - power((STATES(:,1)+45.0000)./12.0000, 2.00000))+0.0400000./(1.00000+exp(( - STATES(:,1)+25.0000)./25.0000))+0.0150000./(1.00000+exp((STATES(:,1)+75.0000)./25.0000))+0.00170000;
    RATES(:,8) = (ALGEBRAIC(:,5) - STATES(:,8))./ALGEBRAIC(:,19);
    ALGEBRAIC(:,6) = 1.00000./(1.00000+exp((STATES(:,1)+26.7000)./5.40000));
    ALGEBRAIC(:,20) =  0.0410000.*exp( - power((STATES(:,1)+47.0000)./12.0000, 2.00000))+0.0800000./(1.00000+exp((STATES(:,1)+55.0000)./ - 5.00000))+0.0150000./(1.00000+exp((STATES(:,1)+75.0000)./25.0000))+0.00170000;
    RATES(:,9) = (ALGEBRAIC(:,6) - STATES(:,9))./ALGEBRAIC(:,20);
    ALGEBRAIC(:,21) = 1.00000./( 45.1600.*exp( 0.0357700.*(STATES(:,1)+50.0000))+ 98.9000.*exp(  - 0.100000.*(STATES(:,1)+38.0000)));
    ALGEBRAIC(:,8) = 1.00000./(1.00000+exp((STATES(:,1)+12.5000)./ - 7.70000));
    RATES(:,12) = (ALGEBRAIC(:,8) - STATES(:,12))./ALGEBRAIC(:,21);
    ALGEBRAIC(:,22) = 3.00000./( 45.1600.*exp( 0.0357700.*(STATES(:,1)+50.0000))+ 98.9000.*exp(  - 0.100000.*(STATES(:,1)+38.0000)));
    ALGEBRAIC(:,10) = 1.00000./(1.00000+exp((STATES(:,1)+12.5000)./ - 7.70000));
    RATES(:,14) = (ALGEBRAIC(:,10) - STATES(:,14))./ALGEBRAIC(:,22);
    ALGEBRAIC(:,23) = 1.00000./( 45.1600.*exp( 0.0357700.*(STATES(:,1)+50.0000))+ 98.9000.*exp(  - 0.100000.*(STATES(:,1)+38.0000)));
    ALGEBRAIC(:,11) = 1.00000./(1.00000+exp((STATES(:,1)+12.5000)./ - 7.70000));
    RATES(:,15) = (ALGEBRAIC(:,11) - STATES(:,15))./ALGEBRAIC(:,23);
    ALGEBRAIC(:,24) = 1.00000./( 0.118850.*exp((STATES(:,1)+80.0000)./28.3700)+ 0.562300.*exp((STATES(:,1)+80.0000)./ - 14.1900));
    ALGEBRAIC(:,13) = 1.00000./(1.00000+exp((STATES(:,1)+138.600)./10.4800));
    RATES(:,17) = (ALGEBRAIC(:,13) - STATES(:,17))./ALGEBRAIC(:,24);
    ALGEBRAIC(:,28) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./CONSTANTS(:,3)).*log(CONSTANTS(:,14)./STATES(:,11));
    ALGEBRAIC(:,29) =  CONSTANTS(:,13).*STATES(:,12).*STATES(:,13).*(STATES(:,1) - ALGEBRAIC(:,28));
    ALGEBRAIC(:,30) =  CONSTANTS(:,15).*STATES(:,14).*(STATES(:,1) - ALGEBRAIC(:,28));
    ALGEBRAIC(:,32) = ( (48.0000./(exp((STATES(:,1)+37.0000)./25.0000)+exp((STATES(:,1)+37.0000)./ - 25.0000))+10.0000).*0.00100000)./(1.00000+exp((STATES(:,1) - (ALGEBRAIC(:,28)+76.7700))./ - 17.0000))+( CONSTANTS(:,17).*(STATES(:,1) - (ALGEBRAIC(:,28)+1.73000)))./( (1.00000+exp(( 1.61300.*CONSTANTS(:,3).*(STATES(:,1) - (ALGEBRAIC(:,28)+1.73000)))./( CONSTANTS(:,1).*CONSTANTS(:,2)))).*(1.00000+exp((CONSTANTS(:,14) - 0.998800)./ - 0.124000)));
    ALGEBRAIC(:,40) = ( (( (( CONSTANTS(:,25).*1.00000)./(1.00000+ 0.124500.*exp((  - 0.100000.*STATES(:,1).*CONSTANTS(:,3))./( CONSTANTS(:,1).*CONSTANTS(:,2)))+ 0.0365000.*CONSTANTS(:,70).*exp((  - STATES(:,1).*CONSTANTS(:,3))./( CONSTANTS(:,1).*CONSTANTS(:,2))))).*CONSTANTS(:,14))./(CONSTANTS(:,14)+CONSTANTS(:,26))).*1.00000)./(1.00000+power(CONSTANTS(:,27)./STATES(:,2), 1.50000));
    ALGEBRAIC(:,31) =  CONSTANTS(:,16).*STATES(:,15).*STATES(:,16).*(STATES(:,1) - ALGEBRAIC(:,28));
    ALGEBRAIC(:,34) =  CONSTANTS(:,18).*STATES(:,17).*CONSTANTS(:,69).*(STATES(:,1) - ALGEBRAIC(:,28));
    ALGEBRAIC(:,38) =  CONSTANTS(:,22).*(STATES(:,1) - ALGEBRAIC(:,28));
    RATES(:,11) = (  - (ALGEBRAIC(:,30)+ALGEBRAIC(:,38)+ALGEBRAIC(:,29)+ALGEBRAIC(:,31)+ALGEBRAIC(:,32)+ALGEBRAIC(:,34)+ ALGEBRAIC(:,40).* - 2.00000).*1.00000)./( CONSTANTS(:,56).*CONSTANTS(:,3));
    ALGEBRAIC(:,25) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./CONSTANTS(:,3)).*log(CONSTANTS(:,9)./STATES(:,2));
    ALGEBRAIC(:,26) =  CONSTANTS(:,8).*power(STATES(:,3), 3.00000).*STATES(:,4).*STATES(:,5).*(STATES(:,1) - ALGEBRAIC(:,25));
    ALGEBRAIC(:,27) =  CONSTANTS(:,10).*STATES(:,7).*( (0.900000+STATES(:,10)./10.0000).*STATES(:,8)+ (0.100000 - STATES(:,10)./10.0000).*STATES(:,9)).*(STATES(:,1) - CONSTANTS(:,11));
    ALGEBRAIC(:,33) =  CONSTANTS(:,18).*STATES(:,17).*CONSTANTS(:,19).*(STATES(:,1) - ALGEBRAIC(:,25));
    ALGEBRAIC(:,35) = ALGEBRAIC(:,33)+ALGEBRAIC(:,34);
    ALGEBRAIC(:,36) =  CONSTANTS(:,20).*(STATES(:,1) - ALGEBRAIC(:,25));
    ALGEBRAIC(:,37) =  CONSTANTS(:,21).*(STATES(:,1) - CONSTANTS(:,23));
    ALGEBRAIC(:,39) = ALGEBRAIC(:,36)+ALGEBRAIC(:,37)+ALGEBRAIC(:,38);
    ALGEBRAIC(:,42) = ( CONSTANTS(:,29).*( power(STATES(:,2), 3.00000).*CONSTANTS(:,24).*exp( 0.0374300.*STATES(:,1).*CONSTANTS(:,31)) -  power(CONSTANTS(:,9), 3.00000).*STATES(:,18).*exp( 0.0374300.*STATES(:,1).*(CONSTANTS(:,31) - 1.00000))))./(1.00000+ CONSTANTS(:,30).*( STATES(:,18).*power(CONSTANTS(:,9), 3.00000)+ CONSTANTS(:,24).*power(STATES(:,2), 3.00000)));
    ALGEBRAIC(:,41) = ( CONSTANTS(:,28).*STATES(:,18))./(STATES(:,18)+0.000400000);
    ALGEBRAIC(:,14) = piecewise({VOI -  floor(VOI./CONSTANTS(:,5)).*CONSTANTS(:,5)>=0.00000&VOI -  floor(VOI./CONSTANTS(:,5)).*CONSTANTS(:,5)<=CONSTANTS(:,6), CONSTANTS(:,7) }, 0.00000);
    RATES(:,1) =  - (ALGEBRAIC(:,26)+ALGEBRAIC(:,27)+ALGEBRAIC(:,29)+ALGEBRAIC(:,30)+ALGEBRAIC(:,31)+ALGEBRAIC(:,35)+ALGEBRAIC(:,32)+ALGEBRAIC(:,39)+ALGEBRAIC(:,40)+ALGEBRAIC(:,42)+ALGEBRAIC(:,41)+ALGEBRAIC(:,14))./CONSTANTS(:,4);
    RATES(:,2) = (  - (ALGEBRAIC(:,26)+ALGEBRAIC(:,36)+ ALGEBRAIC(:,42).*3.00000+ ALGEBRAIC(:,40).*3.00000+ALGEBRAIC(:,33)).*1.00000)./( CONSTANTS(:,56).*CONSTANTS(:,3));
    ALGEBRAIC(:,44) = power(STATES(:,18)./CONSTANTS(:,41), CONSTANTS(:,46));
    ALGEBRAIC(:,45) = power(STATES(:,24)./CONSTANTS(:,42), CONSTANTS(:,47));
    ALGEBRAIC(:,46) = ( CONSTANTS(:,45).*( CONSTANTS(:,43).*ALGEBRAIC(:,44) -  CONSTANTS(:,44).*ALGEBRAIC(:,45)))./(1.00000+ALGEBRAIC(:,44)+ALGEBRAIC(:,45));
    ALGEBRAIC(:,48) = (STATES(:,24) - STATES(:,23))./CONSTANTS(:,48);
    RATES(:,24) = ( ALGEBRAIC(:,46).*CONSTANTS(:,56))./CONSTANTS(:,58) - ( ALGEBRAIC(:,48).*CONSTANTS(:,57))./CONSTANTS(:,58);
    ALGEBRAIC(:,49) =  CONSTANTS(:,52).*STATES(:,18).*(CONSTANTS(:,50) - STATES(:,25)) -  CONSTANTS(:,53).*STATES(:,25);
    RATES(:,25) = ALGEBRAIC(:,49);
    ALGEBRAIC(:,43) =  CONSTANTS(:,32).*(STATES(:,19)+STATES(:,20)).*(STATES(:,23) - STATES(:,6));
    ALGEBRAIC(:,47) = (STATES(:,6) - STATES(:,18))./CONSTANTS(:,49);
    ALGEBRAIC(:,50) = 1.00000./(1.00000+( CONSTANTS(:,63).*CONSTANTS(:,60))./power(CONSTANTS(:,60)+STATES(:,6), 2.00000));
    RATES(:,6) =  ALGEBRAIC(:,50).*((( ALGEBRAIC(:,43).*CONSTANTS(:,57))./CONSTANTS(:,59) - ( ALGEBRAIC(:,47).*CONSTANTS(:,56))./CONSTANTS(:,59)) - ( ALGEBRAIC(:,27).*1.00000)./( 2.00000.*CONSTANTS(:,59).*CONSTANTS(:,3)));
    ALGEBRAIC(:,51) = 1.00000./(1.00000+( CONSTANTS(:,64).*CONSTANTS(:,61))./power(CONSTANTS(:,61)+STATES(:,23), 2.00000));
    RATES(:,23) =  ALGEBRAIC(:,51).*(ALGEBRAIC(:,48) - ALGEBRAIC(:,43));
    ALGEBRAIC(:,52) =  CONSTANTS(:,54).*STATES(:,18).*(CONSTANTS(:,51) - STATES(:,26)) -  CONSTANTS(:,55).*STATES(:,26);
    RATES(:,26) = ALGEBRAIC(:,52);
    ALGEBRAIC(:,53) = ALGEBRAIC(:,49)+ALGEBRAIC(:,52);
    ALGEBRAIC(:,54) = 1.00000./(1.00000+( CONSTANTS(:,63).*CONSTANTS(:,60))./power(CONSTANTS(:,60)+STATES(:,18), 2.00000)+( CONSTANTS(:,65).*CONSTANTS(:,62))./power(CONSTANTS(:,62)+STATES(:,18), 2.00000));
    RATES(:,18) =  ALGEBRAIC(:,54).*(ALGEBRAIC(:,47) - (ALGEBRAIC(:,46)+ALGEBRAIC(:,53)+( ((ALGEBRAIC(:,37) -  2.00000.*ALGEBRAIC(:,42))+ALGEBRAIC(:,41)).*1.00000)./( 2.00000.*CONSTANTS(:,56).*CONSTANTS(:,3))));
   RATES = RATES';
end

% Calculate algebraic variables
function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI)
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        utilOnes = 1;
    else
        statesRowCount = statesSize(1);
        utilOnes = ones(statesRowCount, 1);
    end
    ALGEBRAIC(:,7) = 1.00000./(1.00000+STATES(:,6)./0.0100000);
    ALGEBRAIC(:,9) = 1.00000./(1.00000+exp((STATES(:,1)+24.8000)./3.50000));
    ALGEBRAIC(:,12) = 1.00000./(1.00000+exp((STATES(:,1)+37.6000)./5.90000));
    ALGEBRAIC(:,1) = 1.00000./(1.00000+exp((STATES(:,1)+45.0000)./ - 6.50000));
    ALGEBRAIC(:,15) = 0.00136000./(( 0.320000.*(STATES(:,1)+47.1300))./(1.00000 - exp(  - 0.100000.*(STATES(:,1)+47.1300)))+ 0.0800000.*exp( - STATES(:,1)./11.0000));
    ALGEBRAIC(:,2) = 1.00000./(1.00000+exp((STATES(:,1)+76.1000)./6.07000));
    ALGEBRAIC(:,16) = piecewise({STATES(:,1)>= - 40.0000,  0.000453700.*(1.00000+exp( - (STATES(:,1)+10.6600)./11.1000)) }, 0.00349000./( 0.135000.*exp( - (STATES(:,1)+80.0000)./6.80000)+ 3.56000.*exp( 0.0790000.*STATES(:,1))+ 310000..*exp( 0.350000.*STATES(:,1))));
    ALGEBRAIC(:,3) = 1.00000./(1.00000+exp((STATES(:,1)+76.1000)./6.07000));
    ALGEBRAIC(:,17) = piecewise({STATES(:,1)>= - 40.0000, ( 0.0116300.*(1.00000+exp(  - 0.100000.*(STATES(:,1)+32.0000))))./exp(  - 2.53500e-07.*STATES(:,1)) }, 0.00349000./( ((STATES(:,1)+37.7800)./(1.00000+exp( 0.311000.*(STATES(:,1)+79.2300)))).*(  - 127140..*exp( 0.244400.*STATES(:,1)) -  3.47400e-05.*exp(  - 0.0439100.*STATES(:,1)))+( 0.121200.*exp(  - 0.0105200.*STATES(:,1)))./(1.00000+exp(  - 0.137800.*(STATES(:,1)+40.1400)))));
    ALGEBRAIC(:,4) = 1.00000./(1.00000+exp((STATES(:,1)+15.3000)./ - 5.00000));
    ALGEBRAIC(:,18) =  0.00305000.*exp(  - 0.00450000.*power(STATES(:,1)+7.00000, 2.00000))+ 0.00105000.*exp(  - 0.00200000.*power(STATES(:,1) - 18.0000, 2.00000))+0.000250000;
    ALGEBRAIC(:,5) = 1.00000./(1.00000+exp((STATES(:,1)+26.7000)./5.40000));
    ALGEBRAIC(:,19) =  0.105000.*exp( - power((STATES(:,1)+45.0000)./12.0000, 2.00000))+0.0400000./(1.00000+exp(( - STATES(:,1)+25.0000)./25.0000))+0.0150000./(1.00000+exp((STATES(:,1)+75.0000)./25.0000))+0.00170000;
    ALGEBRAIC(:,6) = 1.00000./(1.00000+exp((STATES(:,1)+26.7000)./5.40000));
    ALGEBRAIC(:,20) =  0.0410000.*exp( - power((STATES(:,1)+47.0000)./12.0000, 2.00000))+0.0800000./(1.00000+exp((STATES(:,1)+55.0000)./ - 5.00000))+0.0150000./(1.00000+exp((STATES(:,1)+75.0000)./25.0000))+0.00170000;
    ALGEBRAIC(:,21) = 1.00000./( 45.1600.*exp( 0.0357700.*(STATES(:,1)+50.0000))+ 98.9000.*exp(  - 0.100000.*(STATES(:,1)+38.0000)));
    ALGEBRAIC(:,8) = 1.00000./(1.00000+exp((STATES(:,1)+12.5000)./ - 7.70000));
    ALGEBRAIC(:,22) = 3.00000./( 45.1600.*exp( 0.0357700.*(STATES(:,1)+50.0000))+ 98.9000.*exp(  - 0.100000.*(STATES(:,1)+38.0000)));
    ALGEBRAIC(:,10) = 1.00000./(1.00000+exp((STATES(:,1)+12.5000)./ - 7.70000));
    ALGEBRAIC(:,23) = 1.00000./( 45.1600.*exp( 0.0357700.*(STATES(:,1)+50.0000))+ 98.9000.*exp(  - 0.100000.*(STATES(:,1)+38.0000)));
    ALGEBRAIC(:,11) = 1.00000./(1.00000+exp((STATES(:,1)+12.5000)./ - 7.70000));
    ALGEBRAIC(:,24) = 1.00000./( 0.118850.*exp((STATES(:,1)+80.0000)./28.3700)+ 0.562300.*exp((STATES(:,1)+80.0000)./ - 14.1900));
    ALGEBRAIC(:,13) = 1.00000./(1.00000+exp((STATES(:,1)+138.600)./10.4800));
    ALGEBRAIC(:,28) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./CONSTANTS(:,3)).*log(CONSTANTS(:,14)./STATES(:,11));
    ALGEBRAIC(:,29) =  CONSTANTS(:,13).*STATES(:,12).*STATES(:,13).*(STATES(:,1) - ALGEBRAIC(:,28));
    ALGEBRAIC(:,30) =  CONSTANTS(:,15).*STATES(:,14).*(STATES(:,1) - ALGEBRAIC(:,28));
    ALGEBRAIC(:,32) = ( (48.0000./(exp((STATES(:,1)+37.0000)./25.0000)+exp((STATES(:,1)+37.0000)./ - 25.0000))+10.0000).*0.00100000)./(1.00000+exp((STATES(:,1) - (ALGEBRAIC(:,28)+76.7700))./ - 17.0000))+( CONSTANTS(:,17).*(STATES(:,1) - (ALGEBRAIC(:,28)+1.73000)))./( (1.00000+exp(( 1.61300.*CONSTANTS(:,3).*(STATES(:,1) - (ALGEBRAIC(:,28)+1.73000)))./( CONSTANTS(:,1).*CONSTANTS(:,2)))).*(1.00000+exp((CONSTANTS(:,14) - 0.998800)./ - 0.124000)));
    ALGEBRAIC(:,40) = ( (( (( CONSTANTS(:,25).*1.00000)./(1.00000+ 0.124500.*exp((  - 0.100000.*STATES(:,1).*CONSTANTS(:,3))./( CONSTANTS(:,1).*CONSTANTS(:,2)))+ 0.0365000.*CONSTANTS(:,70).*exp((  - STATES(:,1).*CONSTANTS(:,3))./( CONSTANTS(:,1).*CONSTANTS(:,2))))).*CONSTANTS(:,14))./(CONSTANTS(:,14)+CONSTANTS(:,26))).*1.00000)./(1.00000+power(CONSTANTS(:,27)./STATES(:,2), 1.50000));
    ALGEBRAIC(:,31) =  CONSTANTS(:,16).*STATES(:,15).*STATES(:,16).*(STATES(:,1) - ALGEBRAIC(:,28));
    ALGEBRAIC(:,34) =  CONSTANTS(:,18).*STATES(:,17).*CONSTANTS(:,69).*(STATES(:,1) - ALGEBRAIC(:,28));
    ALGEBRAIC(:,38) =  CONSTANTS(:,22).*(STATES(:,1) - ALGEBRAIC(:,28));
    ALGEBRAIC(:,25) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./CONSTANTS(:,3)).*log(CONSTANTS(:,9)./STATES(:,2));
    ALGEBRAIC(:,26) =  CONSTANTS(:,8).*power(STATES(:,3), 3.00000).*STATES(:,4).*STATES(:,5).*(STATES(:,1) - ALGEBRAIC(:,25));
    ALGEBRAIC(:,27) =  CONSTANTS(:,10).*STATES(:,7).*( (0.900000+STATES(:,10)./10.0000).*STATES(:,8)+ (0.100000 - STATES(:,10)./10.0000).*STATES(:,9)).*(STATES(:,1) - CONSTANTS(:,11));
    ALGEBRAIC(:,33) =  CONSTANTS(:,18).*STATES(:,17).*CONSTANTS(:,19).*(STATES(:,1) - ALGEBRAIC(:,25));
    ALGEBRAIC(:,35) = ALGEBRAIC(:,33)+ALGEBRAIC(:,34);
    ALGEBRAIC(:,36) =  CONSTANTS(:,20).*(STATES(:,1) - ALGEBRAIC(:,25));
    ALGEBRAIC(:,37) =  CONSTANTS(:,21).*(STATES(:,1) - CONSTANTS(:,23));
    ALGEBRAIC(:,39) = ALGEBRAIC(:,36)+ALGEBRAIC(:,37)+ALGEBRAIC(:,38);
    ALGEBRAIC(:,42) = ( CONSTANTS(:,29).*( power(STATES(:,2), 3.00000).*CONSTANTS(:,24).*exp( 0.0374300.*STATES(:,1).*CONSTANTS(:,31)) -  power(CONSTANTS(:,9), 3.00000).*STATES(:,18).*exp( 0.0374300.*STATES(:,1).*(CONSTANTS(:,31) - 1.00000))))./(1.00000+ CONSTANTS(:,30).*( STATES(:,18).*power(CONSTANTS(:,9), 3.00000)+ CONSTANTS(:,24).*power(STATES(:,2), 3.00000)));
    ALGEBRAIC(:,41) = ( CONSTANTS(:,28).*STATES(:,18))./(STATES(:,18)+0.000400000);
    ALGEBRAIC(:,14) = piecewise({VOI -  floor(VOI./CONSTANTS(:,5)).*CONSTANTS(:,5)>=0.00000&VOI -  floor(VOI./CONSTANTS(:,5)).*CONSTANTS(:,5)<=CONSTANTS(:,6), CONSTANTS(:,7) }, 0.00000);
    ALGEBRAIC(:,44) = power(STATES(:,18)./CONSTANTS(:,41), CONSTANTS(:,46));
    ALGEBRAIC(:,45) = power(STATES(:,24)./CONSTANTS(:,42), CONSTANTS(:,47));
    ALGEBRAIC(:,46) = ( CONSTANTS(:,45).*( CONSTANTS(:,43).*ALGEBRAIC(:,44) -  CONSTANTS(:,44).*ALGEBRAIC(:,45)))./(1.00000+ALGEBRAIC(:,44)+ALGEBRAIC(:,45));
    ALGEBRAIC(:,48) = (STATES(:,24) - STATES(:,23))./CONSTANTS(:,48);
    ALGEBRAIC(:,49) =  CONSTANTS(:,52).*STATES(:,18).*(CONSTANTS(:,50) - STATES(:,25)) -  CONSTANTS(:,53).*STATES(:,25);
    ALGEBRAIC(:,43) =  CONSTANTS(:,32).*(STATES(:,19)+STATES(:,20)).*(STATES(:,23) - STATES(:,6));
    ALGEBRAIC(:,47) = (STATES(:,6) - STATES(:,18))./CONSTANTS(:,49);
    ALGEBRAIC(:,50) = 1.00000./(1.00000+( CONSTANTS(:,63).*CONSTANTS(:,60))./power(CONSTANTS(:,60)+STATES(:,6), 2.00000));
    ALGEBRAIC(:,51) = 1.00000./(1.00000+( CONSTANTS(:,64).*CONSTANTS(:,61))./power(CONSTANTS(:,61)+STATES(:,23), 2.00000));
    ALGEBRAIC(:,52) =  CONSTANTS(:,54).*STATES(:,18).*(CONSTANTS(:,51) - STATES(:,26)) -  CONSTANTS(:,55).*STATES(:,26);
    ALGEBRAIC(:,53) = ALGEBRAIC(:,49)+ALGEBRAIC(:,52);
    ALGEBRAIC(:,54) = 1.00000./(1.00000+( CONSTANTS(:,63).*CONSTANTS(:,60))./power(CONSTANTS(:,60)+STATES(:,18), 2.00000)+( CONSTANTS(:,65).*CONSTANTS(:,62))./power(CONSTANTS(:,62)+STATES(:,18), 2.00000));
end

% Compute result of a piecewise function
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

% Pad out or shorten strings to a set length
function strout = strpad(strin)
    req_length = 160;
    insize = size(strin,2);
    if insize > req_length
        strout = strin(1:req_length);
    else
        strout = [strin, blanks(req_length - insize)];
    end
end
