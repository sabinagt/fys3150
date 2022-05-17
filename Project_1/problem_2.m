%% Plot u(x) n=10
data_10 = readmatrix('exact_solutions_10.csv');
x_10 = data_10(:,1);
u_10 = data_10(:,2);

figure();
%subplot(1,3,1)
scatter(x_10, u_10, 5, 'Filled')
grid on
xlabel('x');
ylabel('u(x)');
title('Plot of u(x) n=10');

%% Plot u(x) n=100
data_100 = readmatrix('exact_solutions_100.csv');
x_100 = data_100(:,1);
u_100 = data_100(:,2);

%subplot(1,3,2)
figure();
scatter(x_100, u_100, 5, 'Filled')
grid on
xlabel('x');
ylabel('u(x)');
title('Plot of u(x) n=100');

%% % Plot u(x) n=1000
data_1000 = readmatrix('exact_solutions_1000.csv');
x_1000 = data_1000(:,1);
u_1000 = data_1000(:,2);

%subplot(1,3,3)
%continuos function
figure();
plot(x_1000, u_1000)
grid on
xlabel('x');
ylabel('u(x)');
title('Plot of u(x)');

