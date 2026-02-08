% Han Zheng
% Beam Analysis Tool

function beam_analysis_tool()
    % units
    fprintf('--- Beam Analysis Tool ---\n');
    units_choice = input('Choose units (1 for N & mm, 2 for lb & inches): ');
    if units_choice == 1
        f_u = 'N'; l_u = 'mm'; m_u = 'N-mm'; s_u = 'MPa';
    else
        f_u = 'lb'; l_u = 'in'; m_u = 'lb-in'; s_u = 'psi';
    end

    % geometry setup
    L = input(['Enter beam length (', l_u, '): ']);
    cs_choice = input('Cross-section (1: Rectangle, 2: Hollow Cylinder): ');
    if cs_choice == 1
        b = input(['Width (', l_u, '): ']);
        h = input(['Height (', l_u, '): ']);
        I = (b * h^3) / 12;
        A = b * h;
        y_max = h / 2;
        Q_max = (b * (h/2) * (h/4)); % First moment of area at neutral axis
        t_min = b;
    else
        Ro = input(['Outer Radius (', l_u, '): ']);
        Ri = input(['Inner Radius (', l_u, '): ']);
        I = (pi/4) * (Ro^4 - Ri^4);
        A = pi * (Ro^2 - Ri^2);
        y_max = Ro;
        Q_max = (2/3) * (Ro^3 - Ri^3); % First moment of area for half-circle
        t_min = 2 * (Ro - Ri);
    end

    % BC
    % 1: Simply Supported (Pinned/Roller), 2: Cantilever
    BC = input('Boundary Condition (1: Simple, 2: Cantilever): ');
    fixed_side = 0; s1_pos = 0; s2_pos = 0;
    if BC == 1
        s1_pos = input(['Position of Left Support (', l_u, '): ']);
        s2_pos = input(['Position of Right Support (', l_u, '): ']);
    else
        fixed_side = input('Fixed End Position (1: Left, 2: Right): ');
    end

    % Loads
    nP = input('Number of Point Forces: ');
    P = zeros(nP, 2);
    for i = 1:nP
        P(i,1) = input(['Force ', num2str(i), ' Mag (', f_u, '): ']);
        P(i,2) = input(['Force ', num2str(i), ' Pos (', l_u, '): ']);
    end

    nM = input('Number of Point Moments: ');
    M_pt = zeros(nM, 2); 
    for i = 1:nM
        M_pt(i,1) = input(['Moment ', num2str(i), ' Mag (', m_u, '): ']);
        M_pt(i,2) = input(['Moment ', num2str(i), ' Pos (', l_u, '): ']);
    end

    nD = input('Number of Distributed Loads: ');
    D = zeros(nD, 3); 
    for i = 1:nD
        D(i,1) = input(['Dist. Load ', num2str(i), ' Mag (', f_u, '/', l_u, '): ']);
        D(i,2) = input(['Dist. Load ', num2str(i), ' Start (', l_u, '): ']);
        D(i,3) = input(['Dist. Load ', num2str(i), ' End (', l_u, '): ']);
    end

    % Calculate Reactions
    R1 = 0; R2 = 0; M_fixed = 0;
    if BC == 1
        sumM_S1 = 0;
        for i=1:nP, sumM_S1 = sumM_S1 + P(i,1)*(P(i,2)-s1_pos); end
        for i=1:nM, sumM_S1 = sumM_S1 + M_pt(i,1); end
        for i=1:nD
            W = D(i,1)*(D(i,3)-D(i,2));
            W_pos = (D(i,3)+D(i,2))/2;
            sumM_S1 = sumM_S1 + W*(W_pos-s1_pos);
        end
        R2 = -sumM_S1 / (s2_pos - s1_pos);
        R1 = -(sum(P(:,1)) + sum(D(:,1).*(D(:,3)-D(:,2))) + R2);
        fprintf('\nReactions: R1 (at %.1f) = %.2f %s, R2 (at %.1f) = %.2f %s\n', s1_pos, R1, f_u, s2_pos, R2, f_u);
    else
        % Cantilever
        R1 = -(sum(P(:,1)) + sum(D(:,1).*(D(:,3)-D(:,2))));
        fixed_pos = 0; if fixed_side == 2, fixed_pos = L; end
        M_sum_fixed = 0;
        for i=1:nP, M_sum_fixed = M_sum_fixed + P(i,1)*(P(i,2)-fixed_pos); end
        for i=1:nM, M_sum_fixed = M_sum_fixed + M_pt(i,1); end
        for i=1:nD
            W = D(i,1)*(D(i,3)-D(i,2));
            W_pos = (D(i,3)+D(i,2))/2;
            M_sum_fixed = M_sum_fixed + W*(W_pos-fixed_pos);
        end
        M_fixed = -M_sum_fixed;
        fprintf('\nReactions: R_fixed = %.2f %s, M_fixed = %.2f %s\n', R1, f_u, M_fixed, m_u);
    end

 
    x = linspace(0, L, 1000);
    V = zeros(size(x));
    M = zeros(size(x));

    for j = 1:length(x)
        curr_x = x(j);
        cv = 0; cm = 0;
        if BC == 1
            if curr_x >= s1_pos, cv = cv + R1; cm = cm + R1*(curr_x - s1_pos); end
            if curr_x >= s2_pos, cv = cv + R2; cm = cm + R2*(curr_x - s2_pos); end
        else
            if fixed_side == 1 
                cv = cv + R1; cm = cm + M_fixed + R1*curr_x;
            end
        end
        for i=1:nP
            if curr_x >= P(i,2)
                cv = cv + P(i,1);
                cm = cm + P(i,1)*(curr_x - P(i,2));
            end
        end
        for i=1:nD
            if curr_x > D(i,2)
                load_end = min(curr_x, D(i,3));
                load_len = load_end - D(i,2);
                W_sub = D(i,1) * load_len;
                W_sub_pos = (load_end + D(i,2))/2;
                cv = cv + W_sub;
                cm = cm + W_sub*(curr_x - W_sub_pos);
            end
        end
        for i=1:nM
            if curr_x >= M_pt(i,2)
                cm = cm + M_pt(i,1);
            end
        end
        V(j) = cv; M(j) = cm;
    end

    % Calculate stress
    sigma = (abs(M) .* y_max) ./ I;
    tau = (abs(V) .* Q_max) ./ (I * t_min);
    
    [max_sig, idx_sig] = max(sigma);
    [max_tau, idx_tau] = max(tau);

    fprintf('\nMax Normal Stress: %.2f %s at x=%.1f %s\n', max_sig, s_u, x(idx_sig), l_u);
    fprintf('Max Shear Stress: %.2f %s at x=%.1f %s\n', max_tau, s_u, x(idx_tau), l_u);

    figure('Color', 'w', 'Name', 'Beam Analysis');
    subplot(2,1,1);
    plot(x, V, 'g', 'LineWidth', 2); grid on;
    ylabel(['Shear (', f_u, ')']); title('Shear Force Diagram');
    line([0 L], [0 0], 'Color', 'k');

    subplot(2,1,2);
    plot(x, M, 'r', 'LineWidth', 2); grid on;
    ylabel(['Moment (', m_u, ')']); xlabel(['Position (', l_u, ')']);
    title('Bending Moment Diagram');
    line([0 L], [0 0], 'Color', 'k');
end