%function [C] = Test_1026(times,sigma,StepSize,k,z)
z=250;Realization=1;StepSize=2*10^(-3);

%{
[C1] = Test_1026(Realization,sqrt(10^(-3)),StepSize,1,z)
[C2] = Test_1026(Realization,sqrt(10^(-3)),StepSize,2,z)
[C3] = Test_1026(Realization,sqrt(10^(-3)),StepSize,4,z)
[C4] = Test_1026(Realization,sqrt(10^(-3)),StepSize,8,z)
[C5] = Test_1025(Realization,sqrt(10^(-3)),StepSize,1,z)
[C6] = Test_1025(Realization,sqrt(10^(-3)),StepSize,2,z)
[C7] = Test_1025(Realization,sqrt(10^(-3)),StepSize,4,z)
[C8] = Test_1025(Realization,sqrt(10^(-3)),StepSize,8,z)
%}
[C5] = Test_1025(Realization,sqrt(10^(-3)),StepSize,1,z)
%[C5] = Test_1029(Realization,sqrt(10^(-3)),StepSize,1,z)




x=1:z;

C9(x)= 3.35;
C10(x)= 3.15;

%plot(x,C10(x),x,C9(x),x,C1(x),x,C2(x),x,C3(x),x,C4(x),x,C5(x),x,C6(x),x,C7(x),x,C8(x))
plot(x,C5(x))
%legend('Wiener co-op','Wiener','C1(co-op)','C2(co-op)','C3(co-op)','C4(co-op)','C5','C6','C7','C8')
xlabel('Training Length') % x-axis label
ylabel('C(bit/channel)') % y-axis label