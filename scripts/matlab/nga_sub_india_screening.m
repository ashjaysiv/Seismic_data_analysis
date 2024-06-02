clc
clear all
filename = 'Copy of NGAsub_MegaFlatfile_RotD50_050_R211022.xlsx';  
[numbers, text, raw] = xlsread(filename);
raw1=string(raw);
raw1(1,:)=[];

% point 1 remove lack of metadata
metadatamag_rem=find(numbers(:,9)==-999);
metadatahyp_rem=find(numbers(:,12)==-999);
metadatajb_rem=find(numbers(:,47)==-999);
metadatavel_rem=find(numbers(:,73)==-999);
metadat_rem=[metadatamag_rem;metadatahyp_rem;metadatajb_rem;metadatavel_rem;];
% metadat_rem=unique(metadat_rem);
raw1(metadat_rem,:)=[];
numbers(metadat_rem,:)=[];


% other than forearc removed
metadata_forearc_flag=find(numbers(:,62)~=2);
% rrup less than 1000 km
metadatarup_rem=find(numbers(:,48)>1000);
% sensor depth
metadatasend_rem=find(numbers(:,71)>2);
% seggregate the interface and intraslab
metadata_inter=find(numbers(:,18)==0 & numbers(:,12)>40);
metadata_intra_slab=find(numbers(:,18)==1 & numbers(:,12)>200);
% T<Tlu
LUT_rem=find(4>numbers(:,42));
%multiple events
metadatamlte_rem=find(numbers(:,20)==1);
% late trigger
metadatalatep_rem=find(numbers(:,40)==1);
% removing PGA>10
metadatapga_rem=find(numbers(:,114)>10);
%% removing ALL POINTS
metadat_rem1=[metadata_forearc_flag;metadatarup_rem;metadata_inter;metadata_intra_slab; LUT_rem;metadatamlte_rem;metadatasend_rem;metadatalatep_rem;metadatapga_rem];
% metadat_rem2=unique(metadat_rem1);
raw1(metadat_rem1,:)=[];
numbers(metadat_rem1,:)=[];

% source review flag
metadatasrc_rev_rem=find(numbers(:,24)==-1); % keep not remove
raw1(metadatasrc_rev_rem,:)=[];
numbers(metadatasrc_rev_rem,:)=[];

% intrface intr slab
metadata_rem_flag1=find(numbers(:,18)==2);
metadata_rem_flag2=find(numbers(:,18)==3);
metadata_rem_flag3=find(numbers(:,18)==4);
metadata_rem_flag4=find(numbers(:,18)==-444);
metadata_rem_flag5=find(numbers(:,18)==-666);
metadata_rem_flag6=find(numbers(:,18)==-777);
metadata_rem_flag7=find(numbers(:,18)==-888);
metadata_rem_flag8=find(numbers(:,18)==-999);
metadat=[metadata_rem_flag1;metadata_rem_flag2;metadata_rem_flag3;metadata_rem_flag4;metadata_rem_flag5;metadata_rem_flag6;metadata_rem_flag7;metadata_rem_flag8];

raw1(metadat,:)=[];
numbers(metadat,:)=[];

data_final=raw1;

% interpolatoion
psa_sub=str2double(data_final(:,117:227));
period = [ 0.01	 0.02	0.03	0.04	0.05	0.075	0.1	 0.12 0.15 0.17 0.2	0.25	0.3	  0.4	 0.5	 0.75	1	2	3	4	];
T=erase(raw(1,117:227),'S');
T=erase(T,'T');
T=str2double(T);
for i=1:20127
psa_sub_final(i,:)=interp1(T,psa_sub(i,:),period);
end

PGA_g=str2double(data_final(:,114));
PGV_cm_sec=str2double(data_final(:,115));
id=find(PGA_g==-999);
PGA_g(id,:)=[];
PGV_cm_sec(id,:)=[];
data_final(id,:)=[];
psa_sub_final(id,:)=[];


% input
target=[log(PGA_g) log(PGV_cm_sec) log(psa_sub_final)];


input(:,1)=str2double(data_final(:,9));
input(:,2)=str2double(data_final(:,47));
input(:,3)=log(str2double(data_final(:,47)));
input(:,4)=log(str2double(data_final(:,73)));
input(:,5)=str2double(data_final(:,12));
input(:,6)=str2double(data_final(:,18));


% eqid
EqID1=str2double(data_final(:,3));
station(:,1)=str2double(data_final(:,69));
station(:,1)=str2double(data_final(:,70));
[p,q,r]=unique(station,'rows');
siteID1=r;


%%indian data
id=find(input_india(:,5)>40);
input_india(id,6)=1;
input_india(~id,6)=0;

input_final=[input; input_india];
target_final=[target;target_india];
eqid=[EqID1;EqID];
siteid=[siteID1;siteID];