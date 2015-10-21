%function [C] = Test_1026(times,sigma,StepSize,k,z)
z=10^(2);Realization=2;StepSize=2*10^(-3);

[C1] = Test_1026(Realization,sqrt(10^(-3)),StepSize,1,z)
[C2] = Test_1026(Realization,sqrt(10^(-3)),StepSize,2,z)
[C3] = Test_1026(Realization,sqrt(10^(-3)),StepSize,4,z)
[C4] = Test_1026(Realization,sqrt(10^(-3)),StepSize,8,z)
[C5] = Test_1025(Realization,sqrt(10^(-3)),StepSize,1,z)
[C6] = Test_1025(Realization,sqrt(10^(-3)),StepSize,2,z)
[C7] = Test_1025(Realization,sqrt(10^(-3)),StepSize,4,z)
[C8] = Test_1025(Realization,sqrt(10^(-3)),StepSize,8,z)

x=1:z;
plot(x,C1(x),x,C2(x),x,C3(x),x,C4(x),x,C5(x),x,C6(x),x,C7(x),x,C8(x))
legend('C1(co-op)','C2(co-op)','C3(co-op)','C4(co-op)','C5','C6','C7','C8')