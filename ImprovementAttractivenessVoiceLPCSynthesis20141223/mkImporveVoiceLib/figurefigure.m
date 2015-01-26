figure(21);
plot(vtlHeightAll_male(:,1),vtlHeightAll_male(:,2),'b.','markersize',12);
grid on;hold on;
plot(vtlHeightAll_female(:,1),vtlHeightAll_female(:,2),'r.','markersize',12);
lsline;
set(gca,'fontsize',24);
xlabel('estimated relative VTL');
ylabel('height (cm)');
title('VTL vs height regression for all ages');
legend('male','female');
saveas(figure(21),'VTL-HeightRegressionForAllAges2.fig');


figure(22);
plot(vtlHeightYoung_male20(:,1),vtlHeightYoung_male20(:,2),'b.','markersize',12);
grid on;hold on;
plot(vtlHeightYoung_female20(:,1),vtlHeightYoung_female20(:,2),'r.','markersize',12);
lsline;
set(gca,'fontsize',24);
xlabel('estimated relative VTL');
ylabel('height (cm)');
title('VTL vs height regression for young');
legend('male','female');
saveas(figure(22),'VTL-HeightRegressionForYoung2.fig');


figure(23);
plot(vtlHeightAdult_male20(:,1),vtlHeightAdult_male20(:,2),'b.','markersize',12);
grid on;hold on;
plot(vtlHeightAdult_female20(:,1),vtlHeightAdult_female20(:,2),'r.','markersize',12);
lsline;
set(gca,'fontsize',24);
xlabel('estimated relative VTL');
ylabel('height (cm)');
title('VTL vs height regression for adults');
legend('male','female');
saveas(figure(23),'VTL-HeightRegressionForAdults2.fig');


figure(24);
plot(vtlWeightAll_male(:,1),vtlWeightAll_male(:,2),'b.','markersize',12);
grid on;hold on;
plot(vtlWeightAll_female(:,1),vtlWeightAll_female(:,2),'r.','markersize',12);
lsline;
set(gca,'fontsize',24);
xlabel('estimated relative VTL');
ylabel('weight (kg)');
title('VTL vs weight regression for all ages');
legend('male','female');
saveas(figure(24),'VTL-WeightRegressionForAllAges2.fig');


figure(25);
plot(vtlWeightYoung_male20(:,1),vtlWeightYoung_male20(:,2),'b.','markersize',12);
grid on;hold on;
plot(vtlWeightYoung_female20(:,1),vtlWeightYoung_female20(:,2),'r.','markersize',12);
lsline;
set(gca,'fontsize',24);
xlabel('estimated relative VTL');
ylabel('weight (kg)');
title('VTL vs weight regression for young');
legend('male','female');
saveas(figure(25),'VTL-WeightRegressionForYoung2.fig');


figure(26);
plot(vtlWeightAdult_male20(:,1),vtlWeightAdult_male20(:,2),'b.','markersize',12);
grid on;hold on;
plot(vtlWeightAdult_female20(:,1),vtlWeightAdult_female20(:,2),'r.','markersize',12);
lsline;
set(gca,'fontsize',24);
xlabel('estimated relative VTL');
ylabel('weight (kg)');
title('VTL vs weight regression for adults');
legend('male','female');
saveas(figure(26),'VTL-WeightRegressionForAdults2.fig');
