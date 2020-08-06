# #############################################################################
# DESCRIPTION:
#
# This shell script converts all .csv file names from Chinese to English
#
# Created by: Peng Wei, peng.wei@geirina.net
# Created on: 01/08/2019
#
# #############################################################################

echo "converting all csv files..."
current_date_time="`date +%Y_%m_%d_%H%M%S`";
path=dump_${current_date_time}
mkdir ${path};
mv *.csv ./${path}

cp ./${path}/LCC换流器表.csv LCC.csv
cp ./${path}/UPFC表.csv  UPFC.csv
cp ./${path}/UPFC串联侧表.csv UPFC_Series.csv
cp ./${path}/UPFC并联侧表.csv UPFC_Parallel.csv
cp ./${path}/VSC换流器表.csv  VSC.csv
#cp ./${path}/刀闸(隔离开关)表.csv  Switch.csv
cp ./${path}/分区表.csv  Zone.csv
cp ./${path}/区域表.csv  Area.csv
cp ./${path}/厂站表.csv  Substation.csv
#cp ./${path}/开关表.csv
cp ./${path}/母线表.csv  Bus.csv
cp ./${path}/节点表.csv  Node.csv
cp ./${path}/负荷表.csv  Load.csv
cp ./${path}/交流线表.csv AC_Line.csv
cp ./${path}/发电机表.csv Generator.csv
# cp ./${path}/数据组表.csv
cp ./${path}/直流线表.csv DC_Line.csv
#cp ./${path}/互感数据表.csv
cp ./${path}/基准容量表.csv  BaseMVA.csv
#cp ./${path}/直流接地表.csv
#cp ./${path}/直流线路表.csv
cp ./${path}/计算方案表.csv  Algorithm.csv
cp ./${path}/马达信息表.csv  Motor.csv
#cp ./${path}/移相变压器表.csv
cp ./${path}/三绕组变压器表.csv Three_winding_transformer.csv
cp ./${path}/两绕组变压器表.csv  Two_winding_transformer.csv
cp ./${path}/计算方案信息表.csv  Algorithm_data.csv
cp ./${path}/串联电容电抗器表.csv Compensator_S.csv
cp ./${path}/并联电容电抗器表.csv Compensator_P.csv
cp ./${path}/静止无功补偿器表.csv SVC.csv

