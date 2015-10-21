%function [C] = Test_1026(times,sigma,StepSize,k,z)
z=10^(2);Realization=2;

[C1] = Test_1026(Realization,sqrt(10^(-3)),10^(-3),1,z)
[C2] = Test_1026(Realization,sqrt(10^(-3)),10^(-3),2,z)
[C3] = Test_1026(Realization,sqrt(10^(-3)),10^(-3),4,z)
[C4] = Test_1026(Realization,sqrt(10^(-3)),10^(-3),8,z)
[C5] = Test_1025(Realization,sqrt(10^(-3)),10^(-3),1,z)
[C6] = Test_1025(Realization,sqrt(10^(-3)),10^(-3),2,z)
[C7] = Test_1025(Realization,sqrt(10^(-3)),10^(-3),4,z)
[C8] = Test_1025(Realization,sqrt(10^(-3)),10^(-3),8,z)

x=1:z;
plot(x,C1(x),x,C2(x),x,C3(x),x,C4(x),x,C5(x),x,C6(x),x,C7(x),x,C8(x))
legend('C1','C2','C3','C4','C5','C6','C7','C8')