%% Initial channel properties 
% included in this section are the physical specifications of the channel
% at the beginning of the simulation and any calculations needed to
% determine essential properties such as mass of elements, flow resistances
% etc...


Qchannel=.150; %MW channel power

H=0; %m height difference from inlet/outlet headers

hmod=1000; %W/m^2.K heat transfer coefficient 

Tenter=100; %C Entering coolant temperature

PSH=210; %kPa supply header pressure

PVH=190; %kPa vent header pressure

Peval=(PSH+PVH)/2/100; % system evaluation pressure

Tmod=60; %C Moderator temperature

Aflow=0.0035; %m^2 Flow area in channel

Dh=0.0074; %m hydraulic diameter of bundle

DPT=0.10338; %m inner pressure tube diameter

DCT=0.12869;%m inner calandria tube thickness

doutclad=0.0138; %m cladding outer diameter

dfuel=0.0122; %m fuel outer diameter

tclad=0.00038; %m thickness of cladding

roc=doutclad/2; %m outer cladding radius

ric=roc-tclad; %m inner cladding radius

Lchannel=5.94; %m length of fuel channel

Lbund=Lchannel; %m length of bundle (for this initial simulation, entire channel will be simulated as one bundle)

tPT=0.00424; %m pressure tube thickness

tCT=0.0014; %m calandria tube thickness

Aipt=pi()*DPT*Lbund;%m^2 inner pressure tube area

Aopt=pi()*(DPT+(2*tPT))*Lbund;%m^2 outer pressure tube area

Aict=pi()*DCT*Lbund;%m^2 inner calandria tube area

Aoct=pi()*(DCT+(2*tCT))*Lbund;%m^2 outer calandria tube area

Afuel=pi()*doutclad*Lbund;%m^2 fuel outer area

Arad=9/37*Afuel; % Area of outer fuel elements that sees PT

doxide=3.00e-6; % thickness of oxide layer (assumed to be low for initial simulation)

sigma=5.670373e-8; %Stefan-Boltzmann constant

eclad=0.325+(0.1246e6*doxide); % dimensionless emissivity value for zirconium cladding
        
ePT=eclad;% emissivity for pressure tube
        
eCT=ePT;% emissivity for calandria tube

ript=DPT/2; %inner pressure tube radius

ropt=ript+tPT; %outer pressure tube radius

rict=DCT/2; % inner calandria tube radius

roct=rict+tCT;% outer calandria tube radius

DhPT=4*Aflow/(pi()*DPT); % pressure tube hydraulic diameter

%% pressure drop calculation based on density difference as well as header pressure

deltaP=(PSH)-PVH+((9.81*H*(XSteam('rhoL_p',Peval)-XSteam('rhoV_p',Peval)))/1000);

%% Resistance Calculation from liquid flow data

DP = [165;
262.5000;
520.5000;
258;
174];  % measured reference pressure drops across the inlet feeder, end fitting, core, and exit end fitting and feeder respectively

M=25.8; % kg/s reference mass flowrate

Rho=780.6; %kg/m^3 coolant density at reference measurements

keff=DP./(M^2); 

reff=keff.*Rho;

RCH=reff(3); %effective resistance of fuel channel

RF=sum(reff(4:5)); %effective resistance of end fitting and feeder




%% mass of zirc in elements
Tref=310; % C reference temperature

TrefK=Tref+273.15;

rhoref=255.66+(0.1024*TrefK); %kg/m^3 reference density

mclad=((doutclad)^2-(doutclad-(2*tclad))^2)/4*Lbund*rhoref; %mass of zirconium in cladding of one element

mPT=((DPT+(2*tPT))^2-(DPT)^2)/4*Lbund*rhoref; %kg mass of zirconium in pressure tube

mCT=((DCT+(2*tCT))^2-(DPT)^2)/4*Lbund*rhoref; %kg mass of zirconium in calandria tube

%% mass of fuel

rhofuel=10970/(1+(2.04e-5*Tref)+(8.7e-9*Tref^2)); %kg/m^3 reference fuel density

mfuel=Lbund*pi()/4*dfuel^2*rhofuel; %kg mass of fuel within one fuel element

%% CO2 thermal conductivity
% as the correlations for CO2 properties are generally very large complex
% and precise mathematical correlations which can take quite a bit of time
% to implement and compute, simple interpolation of measured values will be
% used during this simulation for properties of CO2

% CO2 thermal conductivity for use in interpolation

    kCO2=[0.01051 0.01456 0.01858 0.02257 0.02652 0.03044 0.03814 0.04565 0.05293 0.08491 0.10688 0.11522];

% CO2 thermal conductivity Temperatures for use in interpolation

    kCO2Temp=[-50 0 50 100 150 200 300 400 500 1000 1500 2000];
%% water viscosity above 890C
% as the XSteam package cannot calculate above 900 degrees celsius, these
% va;ues are to continue calculations above vapor temperatures of 900 c

muvap=[43.830
43.868
43.905
43.942
43.979
44.016
44.054
44.091
44.128
44.165
44.202
44.239
44.277
44.314
44.351
44.388
44.425
44.462
44.499
44.536
44.573
44.610
44.647
44.684
44.722
44.759
44.796
44.833
44.869
44.906
44.943
44.980
45.017
45.054
45.091
45.128
45.165
45.202
45.239
45.276
45.312
45.349
45.386
45.423
45.460
45.497
45.533
45.570
45.607
45.644
45.680
45.717
45.754
45.791
45.827
45.864
45.901
45.938
45.974
46.011
46.048
46.084
46.121
46.157
46.194
46.231
46.267
46.304
46.340
46.377
46.414
46.450
46.487
46.523
46.560
46.596
46.633
46.669
46.706
46.742
46.779
46.815
46.852
46.888
46.924
46.961
46.997
47.034
47.070
47.106
47.143
47.179
47.215
47.252
47.288
47.324
47.361
47.397
47.433
47.470
47.506
47.542
47.578
47.615
47.651
47.687
47.723
47.759
47.796
47.832
47.868
47.904
47.940
47.976
48.012
48.049
48.085
48.121
48.157
48.193
48.229
48.265
48.301
48.337
48.373
48.409
48.445
48.481
48.517
48.553
48.589
48.625
48.661
48.697
48.733
48.769
48.805
48.841
48.877
48.912
48.948
48.984
49.020
49.056
49.092
49.128
49.163
49.199
49.235
49.271
49.307
49.342
49.378
49.414
49.450
49.485
49.521
49.557
49.592
49.628
49.664
49.699
49.735
49.771
49.806
49.842
49.878
49.913
49.949
49.984
50.020
50.056
50.091
50.127
50.162
50.198
50.233
50.269
50.304
50.340
50.375
50.411
50.446
50.482
50.517
50.553
50.588
50.623
50.659
50.694
50.730
50.765
50.800
50.836
50.871
50.906
50.942
50.977
51.012
51.048
51.083
51.118
51.154
51.189
51.224
51.259
51.295
51.330
51.365
51.400
51.435
51.471
51.506
51.541
51.576
51.611
51.646
51.681
51.717
51.752
51.787
51.822
51.857
51.892
51.927
51.962
51.997
52.032
52.067
52.102
52.137
52.172
52.207
52.242
52.277
52.312
52.347
52.382
52.417
52.452
52.487
52.522
52.557
52.592
52.626
52.661
52.696
52.731
52.766
52.801
52.835
52.870
52.905
52.940
52.975
53.009
53.044
53.079
53.114
53.148
53.183
53.218
53.253
53.287
53.322
53.357
53.391
53.426
53.461
53.495
53.530
53.565
53.599
53.634
53.668
53.703
53.738
53.772
53.807
53.841
53.876
53.910
53.945
53.979
54.014
54.048
54.083
54.117
54.152
54.186
54.221
54.255
54.289
54.324
54.358
54.393
54.427
54.461
54.496
54.530
54.565
54.599
54.633
54.668
54.702
54.736
54.770
54.805
54.839
54.873
54.908
54.942
54.976
55.010
55.045
55.079
55.113
55.147
55.181
55.216
55.250
55.284
55.318
55.352
55.386
55.420
55.455
55.489
55.523
55.557
55.591
55.625
55.659
55.693
55.727
55.761
55.795
55.829
55.863
55.897
55.931
55.965
55.999
56.033
56.067
56.101
56.135
56.169
56.203
56.237
56.271
56.305
56.338
56.372
56.406
56.440
56.474
56.508
56.541
56.575
56.609
56.643
56.677
56.710
56.744
56.778
56.812
56.846
56.879
56.913
56.947
56.980
57.014
57.048
57.082
57.115
57.149
57.183
57.216
57.250
57.283
57.317
57.351
57.384
57.418
57.451
57.485
57.519
57.552
57.586
57.619
57.653
57.686
57.720
57.753
57.787
57.820
57.854
57.887
57.921
57.954
57.988
58.021
58.054
58.088
58.121
58.155
58.188
58.221
58.255
58.288
58.322
58.355
58.388
58.422
58.455
58.488
58.522
58.555
58.588
58.621
58.655
58.688
58.721
58.754
58.788
58.821
58.854
58.887
58.921
58.954
58.987
59.020
59.053
59.086
59.120
59.153
59.186
59.219
59.252
59.285
59.318
59.351
59.384
59.417
59.451
59.484
59.517
59.550
59.583
59.616
59.649
59.682
59.715
59.748
59.781
59.814
59.847
59.880
59.912
59.945
59.978
60.011
60.044
60.077
60.110
60.143
60.176
60.209
60.241
60.274
60.307
60.340
60.373
60.406
60.438
60.471
60.504
60.537
60.569
60.602
60.635
60.668
60.700
60.733
60.766
60.799
60.831
60.864
60.897
60.929
60.962
60.995
61.027
61.060
61.093
61.125
61.158
61.191
61.223
61.256
61.288
61.321
61.353
61.386
61.419
61.451
61.484
61.516
61.549
61.581
61.614
61.646
61.679
61.711
61.744
61.776
61.809
61.841
61.873
61.906
61.938
61.971
62.003
62.035
62.068
62.100
62.133
62.165
62.197
62.230
62.262
62.294
62.327
62.359
62.391
62.423
62.456
62.488
62.520
62.553
62.585
62.617
62.649
62.682
62.714
62.746
62.778
62.810
62.843
62.875
62.907
62.939
62.971
63.003
63.036
63.068
63.100
63.132
63.164
63.196
63.228
63.260
63.292
63.324
63.356
63.388
63.421
63.453
63.485
63.517
63.549
63.581
63.613
63.645
63.677
63.708
63.740
63.772
63.804
63.836
63.868
63.900
63.932
63.964
63.996
64.028
64.060
64.091
64.123
64.155
64.187
64.219
64.251
64.283
64.314
64.346
64.378
64.410
64.441
64.473
64.505
64.537
64.569
64.600
64.632
64.664
64.695
64.727
64.759
64.791
64.822
64.854
64.886
64.917
64.949
64.981
65.012
65.044
65.075
65.107
65.139
65.170
65.202
65.233
65.265
65.297
65.328
65.360
65.391
65.423
65.454
65.486
65.517
65.549
65.580
65.612
65.643
65.675
65.706
65.738
65.769
65.801
65.832
65.863
65.895
65.926
65.958
65.989
66.020
66.052
66.083
66.115
66.146
66.177
66.209
66.240
66.271
66.303
66.334
66.365
66.396
66.428
66.459
66.490
66.522
66.553
66.584
66.615
66.647
66.678
66.709
66.740
66.771
66.803
66.834
66.865
66.896
66.927
66.958
66.990
67.021
67.052
67.083
67.114
67.145
67.176
67.207
67.238
67.270
67.301
67.332
67.363
67.394
67.425
67.456
67.487
67.518
67.549
67.580
67.611
67.642
67.673
67.704
67.735
67.766
67.797
67.828
67.859
67.889
67.920
67.951
67.982
68.013
68.044
68.075
68.106
68.137
68.167
68.198
68.229
68.260
68.291
68.322
68.352
68.383
68.414
68.445
68.476
68.506
68.537
68.568
68.599
68.629
68.660
68.691
68.722
68.752
68.783
68.814
68.844
68.875
68.906
68.936
68.967
68.998
69.028
69.059
69.090
69.120
69.151
69.182
69.212
69.243
69.273
69.304
69.334
69.365
69.396
69.426
69.457
69.487
69.518
69.548
69.579
69.609
69.640
69.670
69.701
69.731
69.762
69.792
69.823
69.853
69.884
69.914
69.944
69.975
70.005
70.036
70.066
70.096
70.127
70.157
70.188
70.218
70.248
70.279
70.309
70.339
70.370
70.400
70.430
70.461
70.491
70.521
70.551
70.582
70.612
70.642
70.672
70.703
70.733
70.763
70.793
70.824
70.854
70.884
70.914
70.944
70.975]./1000000;

muvaptemp=[890.00
891.00
892.00
893.00
894.00
895.00
896.00
897.00
898.00
899.00
900.00
901.00
902.00
903.00
904.00
905.00
906.00
907.00
908.00
909.00
910.00
911.00
912.00
913.00
914.00
915.00
916.00
917.00
918.00
919.00
920.00
921.00
922.00
923.00
924.00
925.00
926.00
927.00
928.00
929.00
930.00
931.00
932.00
933.00
934.00
935.00
936.00
937.00
938.00
939.00
940.00
941.00
942.00
943.00
944.00
945.00
946.00
947.00
948.00
949.00
950.00
951.00
952.00
953.00
954.00
955.00
956.00
957.00
958.00
959.00
960.00
961.00
962.00
963.00
964.00
965.00
966.00
967.00
968.00
969.00
970.00
971.00
972.00
973.00
974.00
975.00
976.00
977.00
978.00
979.00
980.00
981.00
982.00
983.00
984.00
985.00
986.00
987.00
988.00
989.00
990.00
991.00
992.00
993.00
994.00
995.00
996.00
997.00
998.00
999.00
1000.0
1001.0
1002.0
1003.0
1004.0
1005.0
1006.0
1007.0
1008.0
1009.0
1010.0
1011.0
1012.0
1013.0
1014.0
1015.0
1016.0
1017.0
1018.0
1019.0
1020.0
1021.0
1022.0
1023.0
1024.0
1025.0
1026.0
1027.0
1028.0
1029.0
1030.0
1031.0
1032.0
1033.0
1034.0
1035.0
1036.0
1037.0
1038.0
1039.0
1040.0
1041.0
1042.0
1043.0
1044.0
1045.0
1046.0
1047.0
1048.0
1049.0
1050.0
1051.0
1052.0
1053.0
1054.0
1055.0
1056.0
1057.0
1058.0
1059.0
1060.0
1061.0
1062.0
1063.0
1064.0
1065.0
1066.0
1067.0
1068.0
1069.0
1070.0
1071.0
1072.0
1073.0
1074.0
1075.0
1076.0
1077.0
1078.0
1079.0
1080.0
1081.0
1082.0
1083.0
1084.0
1085.0
1086.0
1087.0
1088.0
1089.0
1090.0
1091.0
1092.0
1093.0
1094.0
1095.0
1096.0
1097.0
1098.0
1099.0
1100.0
1101.0
1102.0
1103.0
1104.0
1105.0
1106.0
1107.0
1108.0
1109.0
1110.0
1111.0
1112.0
1113.0
1114.0
1115.0
1116.0
1117.0
1118.0
1119.0
1120.0
1121.0
1122.0
1123.0
1124.0
1125.0
1126.0
1127.0
1128.0
1129.0
1130.0
1131.0
1132.0
1133.0
1134.0
1135.0
1136.0
1137.0
1138.0
1139.0
1140.0
1141.0
1142.0
1143.0
1144.0
1145.0
1146.0
1147.0
1148.0
1149.0
1150.0
1151.0
1152.0
1153.0
1154.0
1155.0
1156.0
1157.0
1158.0
1159.0
1160.0
1161.0
1162.0
1163.0
1164.0
1165.0
1166.0
1167.0
1168.0
1169.0
1170.0
1171.0
1172.0
1173.0
1174.0
1175.0
1176.0
1177.0
1178.0
1179.0
1180.0
1181.0
1182.0
1183.0
1184.0
1185.0
1186.0
1187.0
1188.0
1189.0
1190.0
1191.0
1192.0
1193.0
1194.0
1195.0
1196.0
1197.0
1198.0
1199.0
1200.0
1201.0
1202.0
1203.0
1204.0
1205.0
1206.0
1207.0
1208.0
1209.0
1210.0
1211.0
1212.0
1213.0
1214.0
1215.0
1216.0
1217.0
1218.0
1219.0
1220.0
1221.0
1222.0
1223.0
1224.0
1225.0
1226.0
1227.0
1228.0
1229.0
1230.0
1231.0
1232.0
1233.0
1234.0
1235.0
1236.0
1237.0
1238.0
1239.0
1240.0
1241.0
1242.0
1243.0
1244.0
1245.0
1246.0
1247.0
1248.0
1249.0
1250.0
1251.0
1252.0
1253.0
1254.0
1255.0
1256.0
1257.0
1258.0
1259.0
1260.0
1261.0
1262.0
1263.0
1264.0
1265.0
1266.0
1267.0
1268.0
1269.0
1270.0
1271.0
1272.0
1273.0
1274.0
1275.0
1276.0
1277.0
1278.0
1279.0
1280.0
1281.0
1282.0
1283.0
1284.0
1285.0
1286.0
1287.0
1288.0
1289.0
1290.0
1291.0
1292.0
1293.0
1294.0
1295.0
1296.0
1297.0
1298.0
1299.0
1300.0
1301.0
1302.0
1303.0
1304.0
1305.0
1306.0
1307.0
1308.0
1309.0
1310.0
1311.0
1312.0
1313.0
1314.0
1315.0
1316.0
1317.0
1318.0
1319.0
1320.0
1321.0
1322.0
1323.0
1324.0
1325.0
1326.0
1327.0
1328.0
1329.0
1330.0
1331.0
1332.0
1333.0
1334.0
1335.0
1336.0
1337.0
1338.0
1339.0
1340.0
1341.0
1342.0
1343.0
1344.0
1345.0
1346.0
1347.0
1348.0
1349.0
1350.0
1351.0
1352.0
1353.0
1354.0
1355.0
1356.0
1357.0
1358.0
1359.0
1360.0
1361.0
1362.0
1363.0
1364.0
1365.0
1366.0
1367.0
1368.0
1369.0
1370.0
1371.0
1372.0
1373.0
1374.0
1375.0
1376.0
1377.0
1378.0
1379.0
1380.0
1381.0
1382.0
1383.0
1384.0
1385.0
1386.0
1387.0
1388.0
1389.0
1390.0
1391.0
1392.0
1393.0
1394.0
1395.0
1396.0
1397.0
1398.0
1399.0
1400.0
1401.0
1402.0
1403.0
1404.0
1405.0
1406.0
1407.0
1408.0
1409.0
1410.0
1411.0
1412.0
1413.0
1414.0
1415.0
1416.0
1417.0
1418.0
1419.0
1420.0
1421.0
1422.0
1423.0
1424.0
1425.0
1426.0
1427.0
1428.0
1429.0
1430.0
1431.0
1432.0
1433.0
1434.0
1435.0
1436.0
1437.0
1438.0
1439.0
1440.0
1441.0
1442.0
1443.0
1444.0
1445.0
1446.0
1447.0
1448.0
1449.0
1450.0
1451.0
1452.0
1453.0
1454.0
1455.0
1456.0
1457.0
1458.0
1459.0
1460.0
1461.0
1462.0
1463.0
1464.0
1465.0
1466.0
1467.0
1468.0
1469.0
1470.0
1471.0
1472.0
1473.0
1474.0
1475.0
1476.0
1477.0
1478.0
1479.0
1480.0
1481.0
1482.0
1483.0
1484.0
1485.0
1486.0
1487.0
1488.0
1489.0
1490.0
1491.0
1492.0
1493.0
1494.0
1495.0
1496.0
1497.0
1498.0
1499.0
1500.0
1501.0
1502.0
1503.0
1504.0
1505.0
1506.0
1507.0
1508.0
1509.0
1510.0
1511.0
1512.0
1513.0
1514.0
1515.0
1516.0
1517.0
1518.0
1519.0
1520.0
1521.0
1522.0
1523.0
1524.0
1525.0
1526.0
1527.0
1528.0
1529.0
1530.0
1531.0
1532.0
1533.0
1534.0
1535.0
1536.0
1537.0
1538.0
1539.0
1540.0
1541.0
1542.0
1543.0
1544.0
1545.0
1546.0
1547.0
1548.0
1549.0
1550.0
1551.0
1552.0
1553.0
1554.0
1555.0
1556.0
1557.0
1558.0
1559.0
1560.0
1561.0
1562.0
1563.0
1564.0
1565.0
1566.0
1567.0
1568.0
1569.0
1570.0
1571.0
1572.0
1573.0
1574.0
1575.0
1576.0
1577.0
1578.0
1579.0
1580.0
1581.0
1582.0
1583.0
1584.0
1585.0
1586.0
1587.0
1588.0
1589.0
1590.0
1591.0
1592.0
1593.0
1594.0
1595.0
1596.0
1597.0
1598.0
1599.0
1600.0
1601.0
1602.0
1603.0
1604.0
1605.0
1606.0
1607.0
1608.0
1609.0
1610.0
1611.0
1612.0
1613.0
1614.0
1615.0
1616.0
1617.0
1618.0
1619.0
1620.0
1621.0
1622.0
1623.0
1624.0
1625.0
1626.0
1627.0
1628.0
1629.0
1630.0
1631.0
1632.0
1633.0
1634.0
1635.0
1636.0
1637.0
1638.0
1639.0
1640.0
1641.0
1642.0
1643.0
1644.0
1645.0
1646.0
1647.0
1648.0
1649.0
1650.0
1651.0
1652.0
1653.0
1654.0
1655.0
1656.0
1657.0
1658.0
1659.0
1660.0
1661.0
1662.0
1663.0
1664.0
1665.0
1666.0
1667.0
1668.0
1669.0
1670.0
1671.0
1672.0
1673.0
1674.0
1675.0
1676.0
1677.0
1678.0
1679.0
1680.0
1681.0
1682.0
1683.0
1684.0
1685.0
1686.0
1687.0
1688.0
1689.0
1690.0
1691.0
1692.0
1693.0
1694.0
1695.0
1696.0
1697.0
1698.0
1699.0
1700.0];

%% Time determination and property matrix creation
% This section sets the length of the simulation and time divisions and
% creates property matrices for the major variables used in the simulation


Time=1000;  %seconds total run time

div=0.1; %s time step length

ind=Time/div; % index to determine time step number

time=zeros(1,(Time/div));

for p=2:length(time)
    
    time(1,p)=time(1,p-1)+div;
end


Tclad=zeros(1,ind); %C cladding temperature



Tvap=zeros(1,ind);%C vapor coolant temperature



TPT=zeros(1,ind);%C pressure tube temperature



Tfuel=zeros(1,ind);%C average fuel temperature



TCT=zeros(1,ind);%C calandria tube temperature



Reynolds=zeros(1,ind);% reynolds number of coolant

Prandtl=zeros(1,ind);% coolant prandtl number

Nusselt=zeros(1,ind);% coolant nusselt number

hcool=zeros(1,ind);% coolant convective heat transfer coefficient 

hrad=zeros(1,ind);%W/m^2.K radiative heat transfer coefficient for cladding/pressure tube

hrpt=zeros(1,ind);% W/m^2.K radiative heat transfer coefficient for calandria tube/pressure tube

kuo2=zeros(1,ind);% W/m.K thermal conductivity of UO2 fuel 

kclad=zeros(1,ind);% W/m.K thermal conductivity of cladding

kPT=zeros(1,ind);% W/m.K pressure tube thermal conductivity

kCT=zeros(1,ind);%W/m.K calandria tube thermal conductivity

kCO2sys=zeros(1,ind);%W/m.K CO2 insulator thermal conductivity

Cpclad=zeros(1,ind);%J/kg.K cladding heat capacity

CpPT=zeros(1,ind);%J/kg.K pressure tube heat capacity

CpCT=zeros(1,ind);%J/kg.K calandria tube heat capacity

Cpvap=zeros(1,ind);%J/kg.K vapor coolant heat capacity

Cpfuel=zeros(1,ind);%J/kg.K fuel heat capacity

hout=zeros(1,ind);% channel exit enthalpy

A1=zeros(1,ind);
    
B1=zeros(1,ind);
    
F1=zeros(1,ind);
    
A2=zeros(1,ind);
    
B2=zeros(1,ind);
    
C2=zeros(1,ind);
    
D2=zeros(1,ind);
       
B3=zeros(1,ind);
    
C3=zeros(1,ind);
    
D3=zeros(1,ind);

F3=zeros(1,ind);

B4=zeros(1,ind);

C4=zeros(1,ind);

D4=zeros(1,ind);

E4=zeros(1,ind);

D5=zeros(1,ind);

E5=zeros(1,ind);

F5=zeros(1,ind);

Rgap=zeros(1,ind); % fuel-cladding gap resistance

R1=zeros(1,ind); 

R2=zeros(1,ind);

R3=zeros(1,ind);

R4=zeros(1,ind);

R5=zeros(1,ind);

CH=zeros(5,ind);

Mcool=zeros(1,ind); %kg mass of coolant within the channel

%% vapor fraction calculation
%This section creates a loop to solve for the vapor fraction (alpha) in the
%channel under low flow conditions

alpha=0;

x=(XSteam('h_pT',Peval,Tenter)-XSteam('hL_p',Peval))/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval));

wsmx=Qchannel*1000/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval));

d=0.01; %dampening factor

delta=0.1;

while delta>=0.0001
    
    wv=(1-alpha)*wsmx/(1+x);
    
    
    hg=XSteam('hV_P',Peval);

    hv=hg+(alpha*Qchannel*1000/wv);
    
    rhov=XSteam('rho_ph',Peval,hv);
    
    rhog=XSteam('rhoV_p',Peval);
    
    rhoave=(rhov+rhog)/2;
    
    a=alpha^2*rhoave/rhog;
    
    alphanew=(1+((1+x)/wsmx*sqrt(deltaP*rhoave/(RCH+(a*RF)))))^-1;
    
    err=alphanew-alpha;
    
    delta=abs(err);
    
    alpha=alpha+(d*err);
    
end
 


%% Initial property calculations
%in this section the initial temperature for the elements and the coolant
%are determined. The temperatures were determined using the assumption that
%saturated fluid was previously flowing in the channel. also calculated is
%the mass flow through the channel

Qbundle=Qchannel; %due to treatment of channel as single bundle to be changed with multiple bundle model

Qel=Qchannel/37; % generation in one element

Qvol=Qel/(pi()/4*dfuel^2*Lchannel); %Volumetric heat generation per pin

mflow=Qbundle*1000*(1-alpha)/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval)); % kg/s 

Tvap(1,1)=XSteam('Tsat_p',Peval); %C initial vapor temperature (saturated)

Reynolds(1,1)=mflow*Dh/Aflow/XSteam('my_pT',Peval,Tvap(1,1)-0.1);% initial reynolds number based on liquid in channel

Prandtl(1,1)=XSteam('Cp_pT',Peval,Tvap(1,1)-0.1)*1000*XSteam('my_pT',Peval,Tvap(1,1)-0.1)/XSteam('tc_pT', Peval,Tvap(1,1)-0.1); % initial prandtl number 

if Reynolds(1,1)<=3000 % initial nusselt number includes consideration for laminar flow
        Nusselt(1,1)=4.36;
else
        
        Nusselt(1,1)=0.023*Reynolds(1,1)^(4/5)*Prandtl(1,1)^0.4;
end

hcool(1,1)=Nusselt(1,1)*XSteam('tcL_p',Peval)/Dh; %J/m^2.K coolant heat transfer coefficient

Tclad(1,1)=(Qel*1000000/Afuel/hcool(1,1))+Tvap(1,1); %C initial cladding temperature 

kuo2(1,1)=kUO2(Tclad(1,1))*1000;%W/m.K initial fuel heat capacity

kclad(1,1)=12.767-(5.4348e-4*(Tclad(1,1)+273.15))+(8.9818e-6*(Tclad(1,1)+273.15)^2); %W/m.K initial cladding heat capacity

Rfuel(1,1)=1/(4*pi()*kuo2(1,1)*Lchannel); % resistance of entire fuel meat

Rclad(1,1)=log(roc/ric)/(2*pi()*kclad(1,1)*Lchannel); %resistance of cladding

Rgap(1,1)=0;% gap resistance to be ignored for now. 

Rconv(1,1)=1/(Afuel/hcool(1,1)); % covective resistance of one element

R1(1,1)=(Rfuel(1,1)/2)+Rgap(1,1)+(Rclad(1,1)/2);% resistance between Tfuel and Tclad

TPT(1,1)=Tvap(1,1); % simplified initial pressure tube temperature

TCT(1,1)=Tmod; % simplified initial calandria tube temperature

Tfuel(1,1)=Tclad(1,1)+(Qel*1000000*R1(1,1)); % initial average fuel temperature

hin=XSteam('hV_p',Peval)*1000; % J/kg enthalpy of steam entering voided section

%% transient calculation
%These calculations determine the five temperatures during the duration of
%the simulation. 

for n=2:ind
    % calculation of vapor heat transfer coefficient. The if loops are to
    % get values if conditions are at saturation conditions 
    
    if Tvap(1,n-1)==XSteam('Tsat_p',Peval)
        
       Reynolds(1,n)=mflow*Dh/Aflow/alpha/XSteam('my_pT',Peval,Tvap(1,1)+0.1);
       
       Prandtl(1,n)=XSteam('Cp_pT',Peval,Tvap(1,1)+0.1)*1000*XSteam('my_pT',Peval,Tvap(1,1)+0.1)/XSteam('tc_pT', Peval,Tvap(1,1)+0.1);
    else
        if Tvap(1,n-1)>890
            
            musys=interp1(muvaptemp,muvap,Tvap(1,n-1));
            
            Reynolds(1,n)=mflow*Dh/Aflow/alpha/musys;
        
            Prandtl(1,n)=XSteam('Cp_pT',Peval,Tvap(1,n-1))*1000*musys/XSteam('tc_pT', Peval,Tvap(1,n-1));
            
            clear musys
            
        else
            Reynolds(1,n)=mflow*Dh/Aflow/alpha/XSteam('my_pT',Peval,Tvap(1,n-1));
        
            Prandtl(1,n)=XSteam('Cp_pT',Peval,Tvap(1,n-1))*1000*XSteam('my_pT',Peval,Tvap(1,n-1))/XSteam('tc_pT', Peval,Tvap(1,n-1));
        end
    end
    
    if Reynolds(1,n)<=3000
        
        Nusselt(1,n)=4.36;
        
    else
        
        Nusselt(1,n)=0.023*(Reynolds(1,n)^(4/5))*(Prandtl(1,n)^0.4);
    end
    
    if Tvap(1,n-1)>XSteam('Tsat_p',Peval)
    
        hcool(1,n)=Nusselt(1,n)*XSteam('tc_pT',Peval,Tvap(1,n-1))/Dh;
    else
        hcool(1,n)=Nusselt(1,n)*XSteam('tc_pT',Peval,Tvap(1,n-1)+0.1)/Dh;
    end
    
    % heat transfer coefficient for radiation- view area of fuel pins is
    % assumed to be equal to the outer half of the 18 outer fuel pins
    hrad(1,n)=sigma*((Tclad(1,n-1)+273.15)+(TPT(1,n-1)+273.15))*(((Tclad(1,n-1)+273.15)^2)+((TPT(1,n-1)+273.15)^2))/((1/eclad)+((1-eclad)/eclad*Arad/Aipt));
    
    hrpt(1,n)=sigma*((TPT(1,n-1)+273.15)+(TCT(1,n-1)+273.15))*(((TPT(1,n-1)+273.15)^2)+((TCT(1,n-1)+273.15)^2))/((1/eclad)+((1-eclad)/eclad*Aopt/Aict));
    
    % heat capacity calculations for main elements
    Cpclad(1,n)=255.66+(0.1024*(Tclad(1,n-1)+273.15)); %J/kg.K
    
    CpPT(1,n)=255.66+(0.1024*(TPT(1,n-1)+273.15)); %J/kg.K
    
    CpCT(1,n)=255.66+(0.1024*(TCT(1,n-1)+273.15));
    
    tau=(Tfuel(1,n-1)+273.15)/1000;
    
    Cpfuel(1,n)=(52.1743+(87.951*tau)-(85.2411*tau^2)+(31.542*tau^3)-(2.6334*tau^4)-(0.71391*tau^-2))/270.03*1000; %J/kg.K
    
    if Tvap(1,n-1)>XSteam('Tsat_p',Peval)
    
        Cpvap(1,n)=XSteam('Cp_pT',Peval,Tvap(1,n-1))*1000;
    else
        Cpvap(1,n)=XSteam('CP_pT',Peval,Tvap(1,n-1)+0.1)*1000;
    end
    % Exit enthalpy calculation based off of temperature of vapor in
    % previous time step
     if Tvap(1,n-1)>XSteam('Tsat_p',Peval)
    
        hout(1,n)=XSteam('h_pT',Peval,Tvap(1,n-1))*1000; % J/kg
    else
        hout(1,n)=XSteam('h_pT',Peval,Tvap(1,n-1)+0.1)*1000;
     end
    
    
     %thermal conductivity of main elements
     
       kuo2(1,n)=kUO2(Tfuel(1,n-1))*1000; % W/m.K
    
    kclad(1,n)=12.767-(5.4348e-4*(Tclad(1,n-1)+273.15))+(8.9818e-6*(Tclad(1,n-1)+273.15)^2); %W/m.K
    
    kPT(1,n)=12.767-(5.4348e-4*(TPT(1,n-1)+273.15))+(8.9818e-6*(TPT(1,n-1)+273.15)^2);
    
    kCT(1,n)=12.767-(5.4348e-4*(TCT(1,n-1)+273.15))+(8.9818e-6*(TCT(1,n-1)+273.15)^2);
    
    kCO2sys(1,n)=interp1(kCO2Temp,kCO2,(Tvap(1,n-1)+Tmod)/2);
    
    % mass of steam in voided section
    
    if Tvap(n-1)==XSteam('Tsat_p',Peval)
        
        Mcool(1,n)=XSteam('rhoV_p',Peval)*Aflow*mflow*div;
    else
        
        Mcool(1,n)=XSteam('rho_pT',Peval,Tvap(1,n-1))*Aflow*mflow*div;
    end
    
    %calculation of resistances. See accompanying document for derivations
    R1(1,n)=(1/8/pi()/kuo2(1,n)/Lchannel)+Rgap(1,n)+(log(roc/ric)/4/pi()/kclad(1,n)/Lchannel);

    R2(1,n)=(log(roc/ric)/4/pi()/kclad(1,n)/Lchannel)+(1/hcool(1,n)/Afuel);
    
    R3(1,n)=(1/hcool(1,n)/Aipt)+(log(ropt/ript)/4/pi()/kPT(1,n)/Lchannel);
    
    R4(1,n)=(log(ropt/ript)/4/pi()/kPT(1,n)/Lchannel)+(log(rict/ropt)/2/pi()/kCO2sys(1,n)/Lchannel)+(log(roct/rict)/4/pi()/kCT(1,n)/Lchannel);
    
    R5(1,n)=(log(roct/rict)/4/pi()/kCT(1,n)/Lchannel)+(1/hmod/Aoct);
    
    %Eq.1 coefficients
    
    A1(1,n)=-1/R1(1,n)/mfuel/Cpfuel(1,n);
    
    B1(1,n)=Tclad(1,n-1)/R1(1,n)/mfuel/Cpfuel(1,n);
    
    F1(1,n)=Qel*1000000/mfuel/Cpfuel(1,n);
    
    %Eq.2 coefficients
    
    A2(1,n)=Tfuel(1,n-1)/mclad/Cpclad(1,n)/R1(1,n);
    
    B2(1,n)=-1/mclad/Cpclad(1,n)*((1/R1(1,n))+(1/R2(1,n))+(hrad(1,n)*Arad));
    
    C2(1,n)=Tvap(1,n-1)/mclad/Cpclad(1,n)/R2(1,n);
    
    D2(1,n)=hrad(1,n)*Arad/mclad/Cpclad(1,n)*TPT(1,n-1);
    
    %Eq.3 coefficients
    
    B3(1,n)=37*Tclad(1,n-1)/Mcool(1,n)/Cpvap(1,n)/R2(1,n);
    
    C3(1,n)=-1/Mcool(1,n)/Cpvap(1,n)*((37/R2(1,n))+(1/R3(1,n)));
    
    D3(1,n)=TPT(1,n-1)/Mcool(1,n)/Cpvap(1,n)/R3(1,n);
    
    F3(1,n)=mflow/alpha/Mcool(1,n)*(hout(1,n)-hin)/Cpvap(1,n);
    
    %Eq.4 coefficients
    
    B4(1,n)=hrad(1,n)*Arad*Tclad(1,n-1)/mPT/CpPT(1,n);
    
    C4(1,n)=Tvap(1,n-1)/mPT/CpPT(1,n)/R3(1,n);
    
    D4(1,n)=-1/mPT/CpPT(1,n)*((1/R3(1,n))+(1/R4(1,n))+(Arad*hrad(1,n))+(hrpt(1,n)*Aopt));
    
    E4(1,n)=TCT(1,n-1)/mPT/CpPT(1,n)*((1/R4(1,n))+(hrpt(1,n)*Aopt));
    
    %Eq.5 coefficients
    
    D5(1,n)=TPT(1,n-1)/mCT/CpCT(1,n)*((1/R4(1,n))+(hrpt(1,n)*Aopt));
    
    E5(1,n)=-1/mCT/CpCT(1,n)*((1/R4(1,n))+(1/R5(1,n))+(hrpt(1,n)*Aopt));
    
    F5(1,n)=Tmod/mCT/CpCT(1,n)/R5(1,n);
    
    % calculation of Temperature values
    
    Tfuel(1,n)=(Tfuel(1,n-1)*exp(A1(1,n)*div))+((1-exp(A1(1,n)*div))*((B1(1,n)+F1(1,n))/-A1(1,n)));
    
    Tclad(1,n)=(Tclad(1,n-1)*exp(B2(1,n)*div))+((1-exp(B2(1,n)*div))*((A2(1,n)+C2(1,n)+D2(1,n))/-B2(1,n)));
    
    Tvap(1,n)=(Tvap(1,n-1)*exp(C3(1,n)*div))+((1-exp(C3(1,n)*div))*((B3(1,n)+D3(1,n)+F3(1,n))/-C3(1,n)));
    
    TPT(1,n)=(TPT(1,n-1)*exp(D4(1,n)*div))+((1-exp(D4(1,n)*div))*((B4(1,n)+C4(1,n)+E4(1,n))/-D4(1,n)));
    
    TCT(1,n)=(TCT(1,n-1)*exp(E5(1,n)*div))+((1-exp(E5(1,n)*div))*((D5(1,n)+F5(1,n))/-E5(1,n)));
    
    if isnan(Tfuel(1,n))
        display('Tfuel')
        break
    elseif isnan(Tclad(1,n))
        display('Tclad')
        break
    elseif isnan(Tvap(1,n))
        display('Tvap')
        break
    elseif isnan(TPT(1,n))
        display('TPT')
        break
    elseif isnan(TCT(1,n))
        display('TCT')
        break
    end
    
end

    Tc=[Tfuel;Tclad;Tvap;TPT;TCT];
    
   plot(time,Tc)
   