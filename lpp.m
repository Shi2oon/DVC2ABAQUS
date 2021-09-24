clear;clc
for iV=1:12
    try
    inO = ['A:\OneDrive - Nexus365\Work\EBSD Data\20-12-03 12 Si\201201_Si3_' ... 
    num2str(iV) '_20kV_15nA_200ms_results\3D_201201_S' num2str(iV) '_Full_map'];
    SN_FE_reinj_dp_xcjin
    clearvars -except iV f
    end
end
f{1} = 'A:\OneDrive - Nexus365\Work\Papers\InSitu Slip\InSitu Slip-bands\3D\S6_1100_10nA_20kV_XEBSD';
f{2} = 'A:\OneDrive - Nexus365\Work\Papers\InSitu Slip\InSitu Slip-bands\3D\S7_1100_10nA_20kV_XEBSD_1';
f{3} = 'A:\OneDrive - Nexus365\Work\Papers\InSitu Slip\InSitu Slip-bands\3D\S7_1100_10nA_20kV_XEBSD_2';
f{4} = 'A:\OneDrive - Nexus365\Work\Papers\InSitu Slip\InSitu Slip-bands\3D\S7_1100_10nA_20kV_XEBSD_3';
for iV=1:4
    try
    inO = f{iV};
    SN_FE_reinj_dp_xcjin
    clearvars -except iV f
    end
end