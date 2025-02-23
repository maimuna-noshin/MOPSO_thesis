function [f, aux] = testFunction(x)
    % Test Function for MOPSO with Battery System and EV Charging
    
    %% Problem-Specific Parameters
    P_load = [200, 180, 150, 120, 100, 90, 80, 100, 150, 200, 250, 300, ...
              350, 400, 450, 500, 480, 450, 400, 350, 300, 250, 200, 180]; % Hourly load (W)
    
    R = [0, 0, 0, 0, 50, 150, 300, 500, 700, 800, 900, 1000, ...
         950, 900, 850, 700, 500, 300, 100, 50, 0, 0, 0, 0]; % Hourly irradiance

    eta_MPPT = 0.9;          % MPPT efficiency
    safety_margin = 0.1;     % Safety margin
    EF = 2.6;                % Emission factor (kg/L)
    dt = 1;                  % Time step (hour)
    T = length(P_load);      % Total time steps
    SOC_min = 0.25;          % Minimum SOC
    SOC_max = 0.85;          % Maximum SOC
    SOC_target = 0.65;       % Target SOC for optimization
    Ec = 0.95;              % Charging efficiency
    Ed = 0.95;              % Discharging efficiency
    battery_unit_capacity = 5; % Battery capacity per unit (kWh)

    % EV Load Demand (0 from 22:00 to 07:00)
    P_EV = [15, 10, 15, 5, 25, 30, 5, 0, 0, 0, 0, 0, 0, 0, 0, 5, 30, 5, 15, 20, 25, 5, 15, 25]; % Hourly EV demand (kW)

    %% Objective Evaluation
    [N, ~] = size(x);
    f = zeros(N, 4);
    
    for i = 1:N
        % Extract and enforce decision variable limits
        N_PV = round(x(i, 1));       % Number of PV panels
        N_Battery = round(x(i, 2));  % Number of battery units
        total_battery_capacity = N_Battery * battery_unit_capacity; % Dynamic capacity

        % Calculate Power Generated
        P_PV = N_PV * (R / 1000) * eta_MPPT * 0.9; % System efficiency included here

        % Initialize Battery SOC and Power
        SOC = zeros(1, T);
        SOC(1) = 0.5;  % Initial SOC (50%)
        P_battery = zeros(1, T);
        LPS = zeros(1, T);
        SOC_dev = 0;   % SOC Deviation Calculation

        % Energy Management and Battery Operation
        for t = 2:T
            P_net = P_PV(t) - P_load(t);
            P_EV_demand = P_EV(t);

            % Battery Control Logic (New SOC Formulas Applied)
            if P_net > 0  % Excess PV power
                if SOC(t-1) < SOC_max
                    P_battery(t) = min(P_net, (SOC_max - SOC(t-1)) * total_battery_capacity / dt / Ec); % Charge
                    SOC(t) = SOC(t-1) + (P_battery(t) * dt * Ec) / total_battery_capacity;
                else
                    P_battery(t) = 0; % Battery Full
                    SOC(t) = SOC(t-1);
                end
            elseif P_net < 0 % Power deficit
                if SOC(t-1) > SOC_min
                    P_battery(t) = max(P_net, (SOC_min - SOC(t-1)) * total_battery_capacity / dt * Ed); % Discharge
                    SOC(t) = SOC(t-1) + (P_battery(t) * dt) / (total_battery_capacity * Ed);
                else
                    P_battery(t) = 0; % Battery Empty
                    SOC(t) = SOC(t-1);
                end
            else
                P_battery(t) = 0;
                SOC(t) = SOC(t-1);
            end

            % 🔹 **Ensure EV gets power if SOC > SOC_min**
            if SOC(t) > SOC_min
                P_EV_supplied = min(P_EV_demand, (SOC(t) - SOC_min) * total_battery_capacity * Ed / dt);
                SOC(t) = SOC(t) - (P_EV_supplied * dt) / (total_battery_capacity * Ed);
            else
                P_EV_supplied = 0;  % No EV charging if below SOC_min
            end

            % 🔹 **Loss of Power Supply Calculation (Only for main load)**
            LPS(t) = max(0, P_load(t) - (P_PV(t) + P_battery(t) - P_EV_supplied));

            % 🔹 **SOC Deviation Calculation (Minimize SOC Fluctuations)**
            SOC_dev = SOC_dev + (SOC(t) - SOC_target)^2;
        end

        % Objective 1: LPSP Calculation
        total_LPS = sum(LPS) * (1 + safety_margin);
        total_load = sum(P_load);
        LPSP = total_LPS / total_load + (0.2 * (1 - SOC_dev));

        % Objective 2: NPC Calculation
        C_PV = 500 * N_PV;       % Example investment cost for PV
        battery_weight = 2; % Increase battery importance
        C_battery = battery_weight * 500 * N_Battery;

        S_V = 0.2 * (C_PV + C_battery);   % Salvage value 
        P_f = 300;         % Penalty factor
        NPC = C_PV + C_battery - S_V + P_f * total_LPS;

        % Objective 3: CO2 Emissions
        fuel_cons = LPS / (0.3 * 1000); % Assuming diesel generator efficiency of 30%
        CO2_emission = sum(fuel_cons) * EF;

        % Assign Objectives (Updated)
        f(i, :) = [NPC, LPSP, CO2_emission, SOC_dev]; % Minimize NPC, LPSP, CO2, and SOC Deviation
    end

    aux = []; % Auxiliary output (empty)
end
