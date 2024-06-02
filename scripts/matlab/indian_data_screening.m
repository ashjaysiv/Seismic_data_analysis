%% indian subduction data
for i=3
    east_comp(1:length(A_EW{i}),i)=A_EW{i};
    north_comp(1:length(A_EW{i}),i)=A_NS{i};
    vert_comp(1:length(A_EW{i}),i)=A_UD{i};
end

for i=3
time(i,1:length(east_comp(:,i)))=0:dt(i,1):(length(east_comp(:,i))-1)*dt(i,1);
end
%% plotting samples
% figure()
% j=120:127;
% for i=1:length(j)
% subplot(3,3,i)
% plot(time(j(i),:),east_comp(:,j(i)),'k')
% end
% 
% figure()
% for i=120:129
% subplot(3,3,i)
% plot(time(i,:),north_comp(:,i),'k')
% end

%% velocity 
for i=1:130
vel_east{i,1}=cumtrapz(time(i,1:length(A_EW{i,:})),A_EW{i,:});
vel_north{i,1}=cumtrapz(time(i,1:length(A_EW{i,:})),A_NS{i,:});
end
%% rotd50
for i=1:130
    theta = 0:1:179;
    east_comp1=A_EW{i,:};
    north_comp1=A_NS{i,:};
    rot_pga = abs(east_comp1*cosd(theta)+north_comp1*sind(theta))';
    pga_rot = max(rot_pga);
    PGA_Rot50(:,i) = median(pga_rot);
    vel_east1=vel_east{i,:}*9.81*100;
    vel_north1=vel_north{i,:}*9.81*100;
    rot_pgv = abs(vel_east1*cosd(theta)+vel_north1*sind(theta))';
    pgv_rot = max(rot_pgv);
    PGV_Rot50(:,i) = median(pgv_rot);
end

%% rotd50 response spectra
for i=1:130
xi=0.05;
period=[0.0100000000000000	0.0200000000000000	0.0300000000000000	0.0400000000000000	0.0500000000000000	0.0750000000000000	0.100000000000000	0.120000000000000	0.150000000000000	0.170000000000000	0.200000000000000	0.250000000000000	0.300000000000000	0.400000000000000	0.500000000000000	0.750000000000000	1	2	3	4];
sPeriod=period;
dt1=dt(i,1);
gacc=A_EW{i,1};
gacc1=A_NS{i,1};
[PSA, PSV, SD, SA, SV, OUT] = responseSpectra(xi, sPeriod, gacc, dt1);
[PSA1, PSV, SD, SA, SV, OUT] = responseSpectra(xi, sPeriod, gacc1, dt1);
psa_north(i,:)=PSA1;
psa_east(i,:)=PSA;
end

for i=88:130
xi=0.05;
period=[0.0100000000000000	0.0200000000000000	0.0300000000000000	0.0400000000000000	0.0500000000000000	0.0750000000000000	0.100000000000000	0.120000000000000	0.150000000000000	0.170000000000000	0.200000000000000	0.250000000000000	0.300000000000000	0.400000000000000	0.500000000000000	0.750000000000000	1	2	3	4];
sPeriod=period;
dt1=dt(i,1);
gacc=A_EW{i,1};
gacc1=A_NS{i,1};
for j=1:length(theta)
temp = gacc*cosd(theta(j))+gacc1*sind(theta(j));
[PSA, PSV, SD, SA, SV, OUT] = responseSpectra(xi, sPeriod, temp, dt1);
RotPSA(j,:) = PSA;
end
PSA_Rot50(i,:) = median(RotPSA);
end

%% inputs

target(:,1)=PGA_Rot50;
target(:,2)=PGV_Rot50;
target(:,3:22)=PSA_Rot50;

input(:,1)=str2double(metadata(:,3)); % mw
input(:,2)=str2double(metadata(:,7)); % rjb
input(:,3)=log(str2double(metadata(:,7))); % log rjb
input(:,5)=str2double(metadata(:,2)); % focal depth

for i=1:130
   if metadata(i,6)=='A'
      input(i,4)=1;
   elseif metadata(i,6)=='B'
       input(i,4)=2;
   elseif metadata(i,6)=='C'
       input(i,4)=3;
   end
end
EqID=str2double(metadata(:,1));

station1(:,1)=metadata(:,4);
station1(:,2)=metadata(:,5);
[p,q,r]=unique(station1,'rows');
siteID=r;


