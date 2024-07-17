%ODE solver for Chronic Stress

[T, Y] = plotOverTime(500, 10, 14);

figure(1)
plot(T,Y,'LineWidth',2);
legend('Glucose','Insulin','\beta cell mass','CRH', 'ACTH', 'GR', 'Cortisol')
xlabel('Hours')
title('Combined Stress \betaIG model')

figure(2)
sgtitle('Combined \beta IG and HPA-axis Model with \rho = 0 in Non-Diabetic Patient')
subplot(2,3,1);
plot(T, Y(:,1), T, Y(:,3),  'LineWidth', 2)
xlabel('Hours')
ylabel('Glucose (mg/dL) and Beta (mg)')
legend('Glucose', '\beta cell mass');
%ylim([70 850]);

subplot(2,3,2);
plot(T, Y(:,2), 'LineWidth', 2)
xlabel('Hours');
ylabel('Insulin (\muU/mL)');
%xlim([0 168]);
%ylim([0 25]);

subplot(2, 3, 3);
plot(T, Y(:, 4), 'LineWidth', 2);
xlabel('Hours');
ylabel('CRH');
% %xlim([0 168]);
% %ylim([0 1.6]);

subplot(2, 3, 4);
plot(T, Y(:, 5), 'LineWidth', 2);
xlabel('Hours');
ylabel('ACTH');
% %xlim([0 168]);
% % ylim([0 0.12]);
% 
subplot(2, 3, 5);
plot(T, Y(:, 6), 'LineWidth', 2);
xlabel('Hours');
ylabel('GR');
% %xlim([0 168]);
% % ylim([0 3]);
% 
subplot(2, 3, 6);
plot(T, Y(:, 7), 'LineWidth', 2);
xlabel('Hours');
ylabel('Cortisol');
% %xlim([0 168]);
% % ylim([0 0.12])

%plotting updated model (trial 1)
figure(3)
subplot(3,3,1)
plot(T, Y(:, 1), 'LineWidth', 2);
xlabel('Hours');
ylabel('Glucose');

subplot(3,3,2)
plot(T, Y(:, 2), 'LineWidth', 2);
xlabel('Hours');
ylabel('Insulin');

subplot(3,3,3)
plot(T, Y(:, 3), 'LineWidth', 2);
xlabel('Hours');
ylabel('Beta');

subplot(3,3,4)
plot(T, Y(:, 4), 'LineWidth', 2);
xlabel('Hours');
ylabel('CRH');

subplot(3,3,5)
plot(T, Y(:, 5), 'LineWidth', 2);
xlabel('Hours');
ylabel('ACTH');

subplot(3,3,6)
plot(T, Y(:, 6), 'LineWidth', 2);
xlabel('Hours');
ylabel('GR');

subplot(3,3,7)
plot(T, Y(:, 7), 'LineWidth', 2);
xlabel('Hours');
ylabel('Cortisol');

subplot(3,3,8)
plot(T, Y(:, 8), 'LineWidth', 2);
xlabel('Hours');
ylabel('\gamma');

subplot(3,3,9)
plot(T, Y(:, 9), 'LineWidth', 2);
xlabel('Hours');
ylabel('\sigma');

% figure(3);
% options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5 1e-5]);
% [~, YSS0] = ode45(@HPA_eq,[0 1000],[1 0.1 0.001 0.1],options, 0, 0);
% c0 = YSS0(end,1); a0 = YSS0(end,2); r0 = YSS0(end,3); o0 = YSS0(end,4);
% 
% options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4 1e-4 1e-4 1e-5 1e-5]);
% [T1,Y1] = ode45(@bigstressfunc,[0 168], [900 0 300 c0 a0 r0 o0], options, 1, 0);
% plot(T1, Y1, 'LineWidth', 2);

function [T, Y] = plotOverTime(total_spikes, duration_of_spike, gap_after_spike)
    f=0; d=0;
    options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5 1e-5]);
    [~, YSS0] = ode45(@HPA_eq,[0 1000],[1 0.1 0.001 0.1],options,f,d);
    c0 = YSS0(end,1); a0 = YSS0(end,2); r0 = YSS0(end,3); o0 = YSS0(end,4);

    T = uint8.empty;
    Y = uint8.empty;

    % INITIAL CONDITIONS
    first_initial = [150 0 800 c0 a0 r0 o0 0 550 0.8];
    % diabtic G_0 = 180
    % non-diabtic G_0 = 100

    for i = 1:total_spikes
        spike_start_time = (duration_of_spike + gap_after_spike) * (i - 1);
        spike_end_time = spike_start_time + duration_of_spike;
        interval_end_time = spike_end_time + gap_after_spike;

        if i ~= 1
            initial = [Y(end, 1) Y(end, 2) Y(end, 3) Y(end, 4) Y(end, 5) Y(end, 6)  Y(end, 7) Y(end, 8)  Y(end, 9) Y(end, 10)];
        else
            initial = first_initial;
        end

        % chronic stress
        f=1; d=0;
        options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4 1e-4 1e-4 1e-5 1e-5 1e-5 1e-5 1e-5]);
        [T1,Y1] = ode45(@bigstressfunc,[spike_start_time spike_end_time], initial, options,f,d);

        second_initial = [Y1(end, 1) Y1(end, 2) Y1(end, 3) Y1(end, 4) Y1(end, 5) Y1(end, 6)  Y1(end, 7) Y1(end, 8)  Y1(end, 9) Y1(end, 10)];
        
        % no stress
        f=0; d=0;
        [T2,Y2] = ode45(@bigstressfunc,[spike_end_time interval_end_time], second_initial, options,f,d);    

        if i == 1
            T = [T1; T2];
            Y = [Y1; Y2];
        else         
            T = [T; T1; T2];
            Y = [Y; Y1; Y2];
        end
    end
end
