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

  cp ./${path}/10型调速器.csv       Governor_10.csv
  cp ./${path}/11-12型调压器.csv    AVR_11_to_12.csv
  cp ./${path}/13型调压器.csv       AVR_13.csv
  cp ./${path}/14型调压器.csv       AVR_14.csv
  cp ./${path}/15型调压器.csv       AVR_15.csv
  cp ./${path}/16型调压器.csv       AVR_16.csv
  cp ./${path}/1型STATCOM.csv      STATCOM_1.csv
  cp ./${path}/1型欠励.csv          UEL_1.csv
  cp ./${path}/1型过励.csv          OEL_1.csv
  cp ./${path}/1型调压器.csv        AVR_1.csv
  cp ./${path}/1型调速器.csv        Governor_1.csv
#  cp ./${path}/1型储能电站.csv
  cp ./${path}/1型光伏发电站.csv     Solar_Generator_1.csv
#  cp ./${path}/1型直流线功率调制.csv
  cp ./${path}/1型静止无功补偿器.csv  SVC_1.csv
  cp ./${path}/2型欠励.csv          UEL_2.csv
  cp ./${path}/2型调压器.csv        AVR_2.csv
  cp ./${path}/2型调速器.csv        Governor_2.csv
  cp ./${path}/2型光伏发电站.csv     Solar_Generator_2.csv
#  cp ./${path}/2型电压源换流器.csv
#  cp ./${path}/2型双馈风力发电机.csv
  cp ./${path}/2型电力系统稳定器.csv  PSS_2.csv
#  cp ./${path}/2型直流线功率调制.csv
#  cp ./${path}/2型直驱风力发电机.csv
  cp ./${path}/2型静止无功补偿器.csv SVC_2.csv
  cp ./${path}/3-10型调压器.csv     AVR_3_to_10.csv
  cp ./${path}/3型调速器.csv        Governor_3.csv
  cp ./${path}/3型电力系统稳定器.csv  PSS_3.csv
#  cp ./${path}/3型直流线功率调制.csv
  cp ./${path}/3型静止无功补偿器.csv  SVC_3.csv
  cp ./${path}/4型调速器.csv         Governor_4.csv
  cp ./${path}/4型电力系统稳定器.csv  PSS_4_6.csv
#  cp ./${path}/4型直流线功率调制.csv
  cp ./${path}/4型静止无功补偿器.csv  SVC_4.csv
  cp ./${path}/5型调速器.csv        Governor_5.csv
  cp ./${path}/5型直流线调节器.csv   DC_Line_Regulator_5.csv
  cp ./${path}/5型电力系统稳定器.csv  PSS_5.csv
  cp ./${path}/5型静止无功补偿器.csv  SVC_5.csv
#  cp ./${path}/6型LCC换流器.csv
  cp ./${path}/6型调速器.csv        Governor_6.csv
  cp ./${path}/7型调速器.csv        Governor_7.csv
  cp ./${path}/7型电力系统稳定器.csv PSS_7.csv
  cp ./${path}/8型调速器.csv        Governor_8.csv
  cp ./${path}/8型电力系统稳定器.csv PSS_8.csv
  cp ./${path}/9型调速器.csv        Governor_9.csv
#  cp ./${path}/Info.csv
#  cp ./${path}/UPFC串联侧模型参数.csv
#  cp ./${path}/UPFC并联侧模型参数.csv
#  cp ./${path}/VSC低频振荡阻尼器.csv
  cp ./${path}/同步机.csv          Synchronous_Machine.csv
#  cp ./${path}/可控串补.csv
#  cp ./${path}/可控高抗.csv
  cp ./${path}/综合负荷.csv      Comprehensive_Load.csv
  cp ./${path}/感应电动机.csv    Induction_Motor.csv
#  cp ./${path}/负荷静特性.csv
#  cp ./${path}/差分方程负荷.csv
  cp ./${path}/直流线调节器.csv    DC_Line_Regulator_1.csv
#  cp ./${path}/配网综合负荷.csv
  cp ./${path}/双馈风力发电机.csv  Doubly_Fed_Wind_Generator.csv
  cp ./${path}/电力系统稳定器.csv  PSS_1.csv
  cp ./${path}/直驱风力发电机.csv  Directly_Driven_Wind_Generator.csv
#  cp ./${path}/柔性直流输电直流线.csv
#  cp ./${path}/调相机附加控制模型.csv
  cp ./${path}/鼠笼异步风力发电机.csv  Squirrel_Cage_Induction_Generator.csv
  cp ./${path}/双馈直驱通用风力发电机.csv Doubly_Fed_Directly_Driven_Wind_Generator.csv
