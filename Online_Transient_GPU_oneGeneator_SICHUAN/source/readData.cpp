/******************************************************************************
 * Copyright (c): (2019-2020) GEIRINA
 * All rights reserved.
 * Project: Power system transient simulation
 *
 * - Author:      Peng Wei, peng.wei@geirina.net
 * - Created on:  Jan. 29, 2020
 * - Last Update: Feb. 15, 2020
 *
 * - This library contains all functions that load data into memory.
 *
*******************************************************************************/

#include "csv.h"
#include "GraphDyn_util.hpp"
#include "readData.hpp"
#include <string>     // std::string, std::to_string

namespace transient_analysis {

void read_bus_data(map<string, BUS> &buses, string path) {
  io::CSVReader<5> in(path);
  in.read_header(io::ignore_extra_column,
                 "name", "v", "ang", "vbase", "island"
                 );
				 /*"Node_Name", "Vm", "Va", "Base_kV" 
				    bus_id name, Vmag, Vang, BasekV*/
  string bus_name;
  real__t Vm, Va, Base_kV; //Vmag //island
  string island;
  real__t system_base = 100.;
  int count = 0; 
  int count_neutral = 0;
  string dm1; //dummy varaibles to skip the Chinese character line
  in.read_row(dm1,dm1,dm1,dm1,dm1);
  while(in.read_row(bus_name, Vm, Va, Base_kV,island)){
    //cout << "DEBUGGING:busname: "<<bus_name << " Vm: " << Vm << endl;
    if(island != "四川.岛0"){//island == "四川.岛-1" 
      //bus_name=="''" //this way also works
      //island=="四川.岛-1" //this way works
      //cout << "SICHUAN DEBUGGING FOUND" << endl;
      continue; // skip buses with island == -1, -1 means dead island, 0 means main island -XZ
    }
    /*if(abs(Base_kV-1.)<EPS){
      ++count_neutral;
      continue;
    }*/
    if( bus_name.find(".1.")!=std::string::npos){
      ++count_neutral; //Disregard neutral point created by topNode.csv(Node.csv)
      continue;
    }

    ++count;
    if (Base_kV < EPS) Base_kV = 1.; //100.0 -WP
    //    if (island != "0") continue;  // skip buses with Valid == -1, -1 means dead island, 0 means main island -XZ
    //if(buses.count(bus_name)>0) continue; //skip same bus
    assert(buses.count(bus_name) == 0);  // make sure bus names are unique

    buses[bus_name].IBase    = system_base / Base_kV;
    buses[bus_name].ZBase    = Base_kV * Base_kV / system_base;
    buses[bus_name].Vm       = Vm;//1.;//Vm;//IN Nodetop.csv, Vm IS PU !!!! Base_kV; //in CIME, Vm is nominal data ----XZ 03/2020
    buses[bus_name].Va       = Va;//0.;//Va; //Va is in degree
    buses[bus_name].Pgen     = 0.;
    buses[bus_name].Qgen     = 0.;
    buses[bus_name].Pload    = 0.;
    buses[bus_name].Qload    = 0.;
    buses[bus_name].Gshunt   = 0.;
    buses[bus_name].Bshunt   = 0.;
    buses[bus_name].invX     = 0.;
    buses[bus_name].bus_type = 1;
    buses[bus_name].bus_name = bus_name; //Matpower
    //slack bus with HARDCODING with bs=3270 for 蜀州1母线
    if (bus_name=="四川.蜀州.500.2649") { //四川.蜀州.500.3270 "3270"  //四川.蜀州.500.2649(20180824年断面)
      cout << "Slack bus found is 四川.蜀州.500.2649" << endl;
      buses[bus_name].bus_type = 3;
    }
  }
  cout << "TopNode Neutral point count is : " << count_neutral << endl;
  cout << "Total buses: " << count << endl;
}

void read_fdpf_data(map<string, BUS> &buses, string path) {  //different from working_demo/eprimodel_sc
  io::CSVReader<5> in(path);
  in.read_header(io::ignore_extra_column,
                 "Bus_Name", "Vm", "Va", "P", "Q"
                 );
  string bus_name;
  real__t Vm, Va, Pn, Qn;
  // int btype;
  
  while(in.read_row(bus_name, Vm, Va, Pn, Qn)) {
    buses[bus_name].Vm = Vm;
    buses[bus_name].Va = Va;
    if (buses[bus_name].bus_type > 1) {
      buses[bus_name].Pgen = -Pn;
      buses[bus_name].Qgen = -Qn;
      buses[bus_name].Pload = 0; //maybe --NO EFFECT
      buses[bus_name].Qload = 0; //maybe --NO EFFECT
    } else {
      buses[bus_name].Pgen = 0; //maybe --NO EFFECT
      buses[bus_name].Qgen = 0; //maybe --NO EFFECT
      buses[bus_name].Pload = Pn;
      buses[bus_name].Qload = Qn;
    }
  }
}
//load.csv
/** loads are PQ buses */
void read_load_data(map<string, BUS> &buses, string path) {
  io::CSVReader<7> in(path);
  in.read_header(io::ignore_extra_column,
                 "node", "off",  "island", "P", "Q", "P_meas", "Q_meas"
                 ); //, "Unit(always divide)"
                 //"Node_Name", "Valid", "Unit", "Pl", "Ql"
  string bus_name;
  int Valid, island;
  real__t Pl, Ql, Pl_given, Ql_given, Pl_meas, Ql_meas;
  int count = 0;
  real__t sumPl=0., sumQl=0.;
  real__t system_base = 100.;//200. //NEW ADDED TO Tigergraph !!! 
  string dm1; //dummy varaibles to skip the Chinese character line
  in.read_row(dm1,dm1,dm1,dm1,dm1, dm1,dm1);
  while(in.read_row(bus_name, Valid, island, Pl_given, Ql_given, Pl_meas, Ql_meas)){
    //if (Valid == 0) continue;  // skip buses with Valid == 0 -from WP
    if (Valid == 1) continue;  // skip buses with Valid == 1, 1 means OFF -XZ
//    if (island == -1) continue;  // skip buses with Valid == -1, -1 means dead island, 0 means main island -XZ
    if (island != 0) continue;  // skip buses with Valid == -1, -1 means dead island, 0 means main island -XZ

    //assert(buses.count(bus_name) > 0); //from wp, if not
    if (buses.count(bus_name) < 1) {
      cout << "Load bus with bs id" << bus_name << "not in the bus bs id list" << endl;
      continue;
    } // skip buses with no match bus bs id
    ++count;
    //if (Unit == 0) {  // unit conversion  -BY XIANG ALWAYS REQUIRED IN EMS DUE TO NO UNIT EXISTED
    /*Pl /= system_base;
    Ql /= system_base; *///direct method
    int PQ_MEASURE_EN=0; //GraphDyn_util definition from EMSModel
    if (PQ_MEASURE_EN==1){ //use measured data
      if(abs(Pl_meas)>EPS) Pl=Pl_meas/system_base;
      else Pl=Pl_given/system_base;
      if(abs(Ql_meas)>EPS) Ql=Ql_meas/system_base;
      else Ql=Ql_given/system_base;      
    }else{  //use given PQ
      Pl=Pl_given/system_base;
      Ql=Ql_given/system_base;
    }

    // if(abs(Pl)>8) Pl*=0.9; //debug purpose before the system can run
    // if(abs(Ql)>5) Ql*=0.9; //debug purpose before the system can run
    //}
    buses[bus_name].bus_type = 1;   // load bus
    buses[bus_name].Pload   += Pl*1.000;
    buses[bus_name].Qload   += Ql*1.000;
    sumPl+=Pl;
    sumQl+=Ql;
  }
  cout << "Total load buses: " << count << ", sumPl=" << sumPl<< " sumQl=" << sumQl << endl;
}

//compensator_P
/** Parallel Capacitors are added to buses as Gshunt and Bshunt */
void read_compensator_P_data(map<string, BUS> &buses, string path) {
  io::CSVReader<5> in(path);
  in.read_header(io::ignore_extra_column,
                 "node", "off", "island", "P", "Q" //Q_rate in TG loaddata
                 );//"Unit(always /base)" Q has to be NEGATIVE
                 //"Node_Name", "Valid", "Unit", "R1", "X1"
  string bus_name;
  int Valid, island;// Unit;
  //real__t R1, X1;
  real__t P, Q, sumcP=0., sumcQ=0.;
  int count=0;
  string dm1; //dummy varaibles to skip the Chinese character line
  in.read_row(dm1,dm1,dm1,dm1,dm1);
  //while(in.read_row(bus_name, Valid, island, R1, X1)){
  while(in.read_row(bus_name, Valid, island, P, Q)){
    ++count;
    //if (Valid == 0) continue;  // skip nonvalid buses
    //if (buses.count(bus_name) == 0) continue;
//    if (abs(R1) > EPS) continue;
    if (Valid == 1) continue;  // skip buses with Valid == 1, 1 means OFF -XZ
//    if (island == -1) continue;  // skip buses with Valid == -1, -1 means dead island, 0 means main island -XZ
    if (island != 0) continue;  // skip buses with Valid == -1, -1 means dead island, 0 means main island -XZ

    //assert(buses.count(bus_name) > 0); //from wp, if not
    if (buses.count(bus_name) < 1) {
      std::cout << "Compensator_P with bs id" << bus_name << "not in the bus bs id list\n";
      continue;
    } // skip buses with no match bus bs id
    P /= 100.;//200.
    Q /= 100.;//200.
    /*if (Unit == 0) {  // unit conversion
      R1 /= buses[bus_name].ZBase;
      X1 /= buses[bus_name].ZBase;
    }
    if (abs(X1) < EPS) X1 = 0.000001;
    real__t denom = R1 * R1 + X1 * X1;
    assert(denom > EPS && abs(X1) > EPS);*/
    int SC_SYS_INPUT=1; //GraphDyn_util definition from EMSModel
    if(SC_SYS_INPUT==0) Q=-Q;
    buses[bus_name].Gshunt += P*1.000; //R1 / denom;
    buses[bus_name].Bshunt += Q*1.000;//-X1 / denom; -Q 4/27 ,it was -Q,but not converge
    //buses[bus_name].invX   += Q; //1. / X1; // NOT USED IN Y_BUS_MATRIX.CPP------------NEEDS REVIEW ---------------- due to compensator_P
    //maybe useful in transient?
    sumcP+=P;
    sumcQ+=Q;

  }
  cout << "Total Compensator count: " << count << endl;
  cout <<  " compP sumP(Gs)=" << sumcP<< " sumQ(Gs)=" << sumcQ << endl;
}

void read_generator_node_data(map<string, BUS> &buses, map<string, GENERATOR> &generators, string path) {
  io::CSVReader<23> in(path);
  in.read_header(io::ignore_extra_column,
                 "off", "island", "node", "P", "Q", "P_meas", "Q_meas",
                 "Ue", "volt", "Gen_Model", "Gen_Par", "AVR_Model",
                 "AVR_Par", "GOV_Model", "GOV_Par", "PSS_Model", "PSS_Par",
                 "Rate_MVA", "Rate_MW", "Xdp", "Xdpp", "X2",
                 "Tj"
                 );/* "Valid", "Mode_Ctrl", "Node_Name", "Pg", "Qg",
                 "V0", "Angle", "Gen_Model", "Gen_Par", "AVR_Model",
                 "AVR_Par", "GOV_Model", "GOV_Par", "PSS_Model", "PSS_Par",
                 "Rate_MVA", "Rate_MW", "Xdp", "Xdpp", "X2",
                 "Tj" */ // P_rate Rate_MVA 
  string bus_name;
  int Valid, island;// bus_type;
  real__t Pg, Qg, V0, Volt, Pg_given, Qg_given, Pg_meas, Qg_meas;
  real__t Gen_Modelx, Gen_Parx, AVR_Modelx, AVR_Parx, GOV_Modelx, GOV_Parx, PSS_Modelx, PSS_Parx;
  int Gen_Model, Gen_Par, AVR_Model, AVR_Par, GOV_Model, GOV_Par, PSS_Model, PSS_Par;
  real__t Rate_MVA, Rate_MW, Xdp, Xdpp, X2, TJ;
  real__t system_base = 100.;//200.
  real__t sumPg=0., sumQg=0., sumPgActive=0., sumPg2Pl=0.;
  int numSlack = 0, numTotal = 0, numValid = 0;
  int count=0;
  string dm1; //dummy varaibles to skip the Chinese character line
  in.read_row(dm1,dm1,dm1,dm1,dm1,  dm1,dm1,dm1,dm1,dm1,  dm1,dm1,dm1,dm1,dm1,  dm1,dm1,dm1,dm1,dm1,  dm1,dm1,dm1);
  while(in.read_row(Valid, island, bus_name, Pg_given, Qg_given, Pg_meas, Qg_meas, V0, Volt,
                    Gen_Modelx, Gen_Parx, AVR_Modelx, AVR_Parx, GOV_Modelx,
                    GOV_Parx, PSS_Modelx, PSS_Parx, Rate_MVA, Rate_MW,
                    Xdp, Xdpp, X2, TJ)){
    /*if (count==0) {
      ++count;
      continue ; //to skip the first Chinese character
    }*/
    Gen_Model = (int)Gen_Modelx; Gen_Par=(int)Gen_Parx; AVR_Model=(int)AVR_Modelx; AVR_Par=(int)AVR_Parx;
    GOV_Model=  (int)GOV_Modelx; GOV_Par=(int)GOV_Parx; PSS_Model=(int)PSS_Modelx; PSS_Par=(int)PSS_Parx;

    ++count;
    //if (Valid == 0) continue;  // skip buses with Valid == 0
//    if (Valid==1 || island == -1) continue; //skip this generator
    if (Valid==1 || island != 0) continue; //skip this generator
    ++numTotal;
    //cout << "Gen busname:" << bus_name << endl;
    if (!(generators.count(bus_name) == 0 && buses.count(bus_name) > 0)){
        cout << "generators.count("<<bus_name<<")="<<generators.count(bus_name) <<endl;
        cout << "buses.count("<<bus_name<<")="<<buses.count(bus_name) <<endl; 
    }
    assert(generators.count(bus_name) == 0 && buses.count(bus_name) > 0);
    /*Pg /= system_base; //Pg is nominal in Unit.csv
    Qg /= system_base; //Qg is nominal in Unit.csv*/
    int PQ_MEASURE_EN=0; //GraphDyn_util definition from EMSModel
    if (PQ_MEASURE_EN==1){ //use measured data
      if(abs(Pg_meas)>EPS) Pg=Pg_meas/system_base;
      else Pg=Pg_given/system_base;
      if(abs(Qg_meas)>EPS) Qg=Qg_meas/system_base;
      else Qg=Qg_given/system_base;    
    }else{  //use given PQ
      // if (Rate_MVA<300.) Pg_given*=1.1;
      Pg=Pg_given/system_base;
      Qg=Qg_given/system_base;
    }
    sumPg+=Pg;
    sumQg+=Qg;

    if (Gen_Model >= 7) {  // all generators that are not implemented are treated as Loads >=37 ---there are 41 generators!! not consider as gentype
      buses[bus_name].bus_type  = 1;   // load
      buses[bus_name].Pload    += -Pg*1.000;
      buses[bus_name].Qload    += -Qg*1.000;
      buses[bus_name].Pgen      = 0.;
      buses[bus_name].Qgen      = 0.;
      sumPg2Pl += Pg;
    } else {    ///////BUS TYPE PROBLEM COMPARISONS!!!!!!!!!!!!!!!!!!!!!!!!!

      /*
      if (bus_type == 0 && numSlack >= 0) {// generator, treated as PV bus (2) or slack bus (3)
        buses[bus_name].bus_type = 3;
        ++numSlack;
        buses[bus_name].Vm = V0; //if bs==3270 !!!!!!!!!!!!!!!!!!!hardcoding
        buses[bus_name].Va = Va;
      } else {
        buses[bus_name].bus_type = 2;
        buses[bus_name].Vm    = V0;
        buses[bus_name].Pgen += Pg;
      }*/ //WEIPENG's METHOD HAS TO CHANGE DUE TO LACK OF bus_type in unit.csv

      //if (bs == '3270' && numSlack >= 0) {// generator, treated as PV bus (2) or slack bus (3)
      //  buses[bus_name].bus_type = 3;
      //  ++numSlack;
      //  buses[bus_name].Vm = V0; //if bs==3270 !!!!!!!!!!!!!!!!!!!HARDCODING --XZ, 03/2020
      //  buses[bus_name].Va = Va;
      //} else {

      // if (Volt <EPS){  //1.EPS, actually useless in Sichuan system since no Volt<EPS
      //   buses[bus_name].bus_type = 2;
      //   buses[bus_name].Vm    = 1.;
      //   buses[bus_name].Pgen += 0.;
      // remarks:|| AVR_Par==196 || AVR_Par==198
      // ||  (Rate_MVA>=12&&Rate_MVA<20)  ||  (Rate_MVA>=0&&Rate_MVA<600) (Rate_MVA>=16&&Rate_MVA<20.)||  
      if (Volt <1. || Pg<EPS  ){ //This is espcially for the case of 四川德昌风电厂。|| AVR_Model!=12|| Pg>6.
        // buses[bus_name].bus_type = 2; //in tigergraph it is type2
        // buses[bus_name].Vm    = 1.;
        // buses[bus_name].Pgen += 0.;
        buses[bus_name].bus_type  = 1;   // load
        buses[bus_name].Pload    += -Pg*1.000;
        buses[bus_name].Qload    += -Qg*1.000;
        buses[bus_name].Pgen      = 0.;
        buses[bus_name].Qgen      = 0.;
        sumPg2Pl += Pg;
      } else {
        buses[bus_name].bus_type = 2;
        buses[bus_name].Vm    = V0/Volt; //V0 is nominal voltage, Volt is base_kV
        buses[bus_name].Pgen += Pg*1.000;
        buses[bus_name].Qgen += Qg*1.000; //no influence on the system power flow
        // buses[bus_name].Va = Va; //do not need due to same Va is included in Node.csv
        sumPgActive+=Pg;

      //}
      ++numValid;
      
      generators[bus_name].bus_name  = bus_name;
      generators[bus_name].Gen_Model = abs(Gen_Model);
      generators[bus_name].Gen_Par   = Gen_Par;
      generators[bus_name].AVR_Model = AVR_Model > 12 ? 1 : AVR_Model;
      generators[bus_name].AVR_Par   = AVR_Par;
      generators[bus_name].GOV_Model = GOV_Model;
      generators[bus_name].GOV_Par   = GOV_Par;
      generators[bus_name].PSS_Model = PSS_Model;//1;//
      generators[bus_name].PSS_Par   = PSS_Par;//2;//
      generators[bus_name].Rate_MVA  = Rate_MVA;
      generators[bus_name].Rate_MW   = Rate_MW;
      generators[bus_name].Xdp       = Xdp < EPS  ? 0.0001 : Xdp;
      generators[bus_name].Xdpp      = Xdpp < EPS ? 0.0001 : Xdpp;
      generators[bus_name].X2        = X2;
      generators[bus_name].TJ        = TJ;

      // extra code for no AVR/PSS/GOV  --does not work well, with more fluctuations
      // if(Gen_Par==1252){ //bus_name == "四川.长河坝厂.20.16957" 
      // // generators[bus_name].Gen_Par   = 40;//81;
      // generators[bus_name].AVR_Model = 12;//12;
      // generators[bus_name].AVR_Par   = 62;//138;//
      // generators[bus_name].GOV_Model = 8;//8;
      // generators[bus_name].GOV_Par   = 25;//59;//
      // generators[bus_name].PSS_Model = 4;//4;
      // generators[bus_name].PSS_Par   = 29;//63;//
      // }
      // if(Gen_Par==40){
      // // generators[bus_name].Gen_Par   = 40;//81;
      // // generators[bus_name].AVR_Model = 12;//12;
      // // generators[bus_name].AVR_Par   = 62;//138;//
      // // generators[bus_name].GOV_Model = 8;//8;
      // // generators[bus_name].GOV_Par   = 25;//59;//
      // generators[bus_name].PSS_Model = 8;//4;
      // generators[bus_name].PSS_Par   = 1;//63;//
      // }
      if( (Rate_MVA>=16&&Rate_MVA<20.) ){
        // generators[bus_name].Gen_Model = 1;
        // generators[bus_name].Gen_Par   = 2;//81;
        //generators[bus_name].AVR_Model = 0;//12;
        //generators[bus_name].AVR_Par   = 0;//138;//
        generators[bus_name].GOV_Model = 8;//8;
        generators[bus_name].GOV_Par   = 59;//59;//
        //generators[bus_name].PSS_Model = 0;//4;
        //generators[bus_name].PSS_Par   = 0;//63;//             
      }   


      }
    }   ///////BUS TYPE PROBLEM COMPARISONS!!!!!!!!!!!!!!!!!!!!!!!!!
  }
  cout << "Total GEN count: " << numTotal << endl;
  cout << "Total valid GEN count: " << numValid << endl;
  cout <<  " sumPg=" << sumPg<< " sumQg=" << sumQg << endl;
  cout <<  " sumPgActive=" << sumPgActive<< " sumPg2Pl=" << sumPg2Pl << endl;
}

void update_line(unordered_map<string, LINE> &lines,
                 string from, string to, real__t R, real__t X,
                 real__t Bh, real__t Gm, real__t Bm, real__t tap, string id_name) {
  if (abs(R) <= EPS && abs(X) <= EPS) X = 1.e-3; //X=1.e-6 -> 1.e-3
  real__t denom = R * R + X * X;
  assert(abs(denom)>EPS);
  string id = from + "#" + to;
  
  if (lines.count(id) == 0) {  // initialize
    lines[id].G      = 0;
    lines[id].B      = 0;
    lines[id].invX   = 0;
    lines[id].B_half = 0;
    lines[id].Gm     = 0;
    lines[id].Bm     = 0;
    lines[id].tap    = 0;
    lines[id].tap_count = 0;
  }
  
  lines[id].from    = from;
  lines[id].to      = to;
  lines[id].G      += +R / denom; //!
  lines[id].B      += -X / denom; //!
  lines[id].invX   += 1. / X;  //!
  lines[id].B_half += Bh;
  lines[id].Gm     += Gm; //!
  lines[id].Bm     += Bm; //!
  lines[id].tap    += tap; //!
  // lines[id].id_name = id_name;
  /*if(lines[id].tap <EPS){
    lines[id].tap = 1.;
  }////////*/
  lines[id].tap_count++;
  //assert( abs(lines[id].G)>EPS);
  assert( abs(lines[id].B)>EPS);
  assert( abs(lines[id].invX)>EPS);
  //if(X<0) cout << " caution line from" << from << " to " << to << "X is NEGATIVE =" << X << endl;
  /*if (abs(lines[id].tap)<EPS){
    cout << "line from" << from << " to " << to << "total tap is " << lines[id].tap << " local tap is " << tap << endl;
  }
  assert( abs(lines[id].tap)>EPS);
//abs(lines[id].Gm)>EPS && abs(lines[id].Bm)>EPS && */
}

/** DC lines are converted to AC lines */
void read_DC_line_data(map<string, BUS> &buses, unordered_map<string, LINE> &lines, string path) {
  io::CSVReader<6> in(path);
  in.read_header(io::ignore_extra_column,
                 "ID_Name", "Valid", "I_Name", "J_Name", "Rl_Ohm",
                 "Ll_mH"
                 );
  string id, from, to;
  int Valid;
  real__t R, L;
  while(in.read_row(id, Valid, from, to, R, L)){
    if (Valid == 0) continue; // lines not connected
    
    R /= buses[from].ZBase;
    L *= 2. * PI * 50. * 0.001 / buses[from].ZBase;
    
    update_line(lines, from, to, R, L, 0., 0., 0., 1.,from); //Matpower, from and to are placeholder, need to change
    update_line(lines, to, from, R, L, 0., 0., 0., -1.,to); //Matpower, from and to are placeholder, need to change
  }
}

//AC_line.csv
void read_AC_line_data(map<string, BUS> &buses, unordered_map<string, LINE> &lines, string path) {
  io::CSVReader<11> in(path);
  in.read_header(io::ignore_extra_column,
                 "id", "I_node", "J_node", "I_off","J_off", 
                 "I_island", "J_island" , "R*", "X*", "B*", "name"
                 );
                 /*"ID_Name", "Valid", "I_Name", "J_Name", "I_Break",
                 "J_Break", "Unit", "R1", "X1", "B1_Half" */
  string id, from, to, id_name;
  int I_Break, J_Break, I_island, J_island;
  real__t R, X, B_half; //sichuan data is B_half already so NO NEED for Bhalf
  int count = 0;
  string dm1; //dummy varaibles to skip the Chinese character line
  in.read_row(dm1,dm1,dm1,dm1,dm1,  dm1,dm1,dm1,dm1,dm1, dm1);
  while(in.read_row(id, from, to, I_Break, J_Break,
                    I_island, J_island, R, X, B_half, id_name)){
//    if (Valid == 0 || I_Break == 0 || J_Break == 0) continue; // line not valid
//    if (I_island == -1 || J_island == -1 || I_Break == 1 || J_Break == 1) continue; // line is neither on (main) island nor online
    if (I_island != 0 || J_island != 0 || I_Break == 1 || J_Break == 1) continue; // line is neither on (main) island nor online
    if (buses.count(from) == 0 || buses.count(to) == 0) {
      cout << "AC bus line from" << from << " to "<< to << "not in the list" << endl;
      continue;
    }
    
    ++count;
    //if (Unit == 0) {    // convert to p.u.
      ////////////////////////////IF R* X* B* means per unitized values //////////////
      /*R /= buses[from].ZBase;
      X /= buses[from].ZBase;
      B_half *= buses[from].ZBase;*/
      ////////////////////////////////////////////////////////////////////////////////
    //}

    if( id_name.find("补装置")!=std::string::npos){
      cout<< "线补测试-- id_name:" << id_name << "x=" << X << "." << endl;
    }

        if(R<=0) {R=-R;}//else{R=0.25*X;} //NEW ADD
    assert(abs(X)>EPS);
    int MFILEPRINT_EN=0; //GraphDyn_util definition from EMSModel
    if (MFILEPRINT_EN==0){
      update_line(lines, from, to, R, X, B_half, 0., 0., 1., id_name);//tap was 1. //Matpower 1. -> 100861. in DEBUGER
      update_line(lines, to, from, R, X, B_half, 0., 0., -1., id_name);//tap was -1. //Matpower//
    }else{
      update_line(lines, from, to, R, X, B_half, 0., 0., 100861., id_name);//tap was 1. //Matpower 1. -> 100861. in DEBUGER
    }
  }
  cout << "AC line count: " << count << endl;
}

void read_two_winding_transformer_data(map<string, BUS> &buses, unordered_map<string, LINE> &lines, string path) {
  io::CSVReader<15> in(path);
  in.read_header(io::ignore_extra_column,
                 "id", "I_island", "J_island", "I_off", "J_off", 
                 "I_node", "J_node", "Ri*", "Xi*", "Rj*", "Xj*","G", "B", "I_t","name"
                 );/*"ID_Name", "Valid", "I_Name", "J_Name", "R1",
                 "X1", "Gm", "Bm", "Tk" */
  string id, from, to, id_name;
  int I_island, J_island, I_off, J_off;
  real__t Ri,Xi,Rj,Xj;
  real__t R, X, Gm, Bm, tap;
  string Gmstr, Bmstr;
  int count = 0;
  while (in.read_row(id, I_island, J_island, I_off, J_off, from, to, Ri, Xi, Rj,Xj,Gmstr, Bmstr, tap, id_name)) {
    /*if (count==-1) { //DO NOT NEED FOR THIS FUNCTION
      ++count;
      continue ; //to skip the first Chinese character
    }*/
//    if (I_island == -1 || J_island == -1 || I_off == 1 || J_off == 1) continue; //skip buses with NOT on main island nor ONLINE
    if (I_island != 0 || J_island != 0 || I_off == 1 || J_off == 1) continue; //skip buses with NOT on main island nor ONLINE

    //if (Valid == 0) continue;  // skip buses with Valid == 0
    if (buses.count(from) == 0 || buses.count(to) == 0) {
      cout << "Two-windingTX from" << from << " to "<< to << "not in the list" << endl;
      continue;}
    ++count;
    if(Gmstr=="''"){
       Gm = 0.;
       //cout << "Gmstr" <<endl;
    }else{
       Gm = stod(Gmstr);
    }
    if(Bmstr=="\'\'"){
       Bm = 0.;
       //cout << "BBBmstr" <<endl;
    }else{
       Bm = stod(Bmstr);
    }
    R = Ri;//+Rj;
    X = Xi;//+Xj;
    if(R<=0) {R=-R;}//else{R=0.25*X;} //NEW ADD
    if(X<=0) {
      X=-X;
      cout << "2TX bus name "<< id_name <<" number" << from << " to " << to << "X is NEGATIVE (fixed to x_new=-x." << endl;
    }//NEW ADD
    if(abs(X)<EPS) {
      cout << "2TX bus number" << from << " to " << to << "X is ZERO" << endl;
      X=0.001;}
      int MFILEPRINT_EN=0; //GraphDyn_util definition from EMSModel
      int SC_SYS_INPUT=1; //GraphDyn_util definition from EMSModel
    if (MFILEPRINT_EN==0){
      if(SC_SYS_INPUT==1) tap=-tap;
      update_line(lines, from, to, R, X, 0., Gm, Bm, tap, id_name);
      update_line(lines, to, from, R, X, 0., Gm, Bm, -tap, id_name);//Matpower//
    }else {
      if(SC_SYS_INPUT==1) tap=-tap; //added on 6/23/2020 to 以前没有-号
      update_line(lines, from, to, R, X, 0., Gm, Bm, tap, id_name);
    }
  }
  cout << "Two winding transformer count: " << count << endl;
}

void read_three_winding_transformer_data(map<string, BUS> &buses, unordered_map<string, LINE> &lines, string path) {
  io::CSVReader<22> in(path);
  in.read_header(io::ignore_extra_column,
                 "id",  "I_island", "K_island", "J_island", "I_off", "K_off", "J_off", 
                 "I_node", "K_node", "J_node",
                 "Ri*", "Xi*", "Rk*", "Xk*",
                 "Rj*", "Xj*", "G", "B", "I_t",
                 "K_t", "J_t", "name"
                 );/* "ID_Name", "Valid", "Name_1", "Name_2", "Name_3",
                 "Name_N", "R1", "X1", "R2", "X2",
                 "R3", "X3", "Gm", "Bm", "Tk1",
                 "Tk2", "Tk3"*/
  //#include <string>     // std::string, std::to_string
  string name1, name2, name3, id_name;// nameN;
  int id;
  //int Valid;
  int I_island, K_island, J_island, I_off, K_off, J_off;
  int count = 0;
  real__t R1, X1, R2, X2, R3, X3, Gm, Bm, tap1, tap2, tap3;
  string Gmstr,Bmstr;
  while (in.read_row(id, I_island, K_island,J_island, I_off, K_off, J_off, name1, name2, name3,
                     R1, X1, R2, X2,
                     R3, X3, Gmstr, Bmstr, tap1,
                     tap2, tap3, id_name)) {
    /*if (count==-1) { //DO NOT NEED FOR THIS FUNCTION
      ++count;
      continue ; //to skip the first Chinese character
    }*/
    //if (I_island == -1 || K_island == -1 || J_island == -1 || I_off == 1 || K_off == 1|| J_off == 1) continue; //skip buses with NOT on main island nor ONLINE
//    if (I_island != 0 || K_island != 0 || J_island != 0 || I_off == 1 || K_off == 1|| J_off == 1) continue; //skip buses with NOT on main island nor ONLINE
    if ((I_off+K_off+J_off)>=2) continue;
		
	  if( !((I_off==0 && I_island == 0)||(J_off==0 && J_island == 0) ||(K_off==0 && K_island == 0)) )continue;
		
    //if (Valid == 0) continue;             // skip when Valid == 0
    std::string nameN = std::to_string((-1)*id);
      if(id==5751){
        cout << ">>> id is |" << nameN <<"|" <<endl;
      }
      if(nameN=="-5751"){
        cout << ">>> '-5751 neutral Ioff='" << I_off << " I_island = " << I_island << endl;
      }  
    // ****************** 3winding Transformer Neutral Point Construction ********
    real__t system_base =100.;
    real__t Base_kV = 1.;
    //a name update to the circuit
    nameN=id_name+".中性点";
    buses[nameN].IBase    = system_base / Base_kV;
    buses[nameN].ZBase    = Base_kV * Base_kV / system_base;
    buses[nameN].Vm       = 1.;//IN Nodetop.csv, Vm IS PU !!!! Base_kV; //in CIME, Vm is nominal data ----XZ 03/2020
    buses[nameN].Va       = 0.; //Va is in degree
    buses[nameN].Pgen     = 0.;
    buses[nameN].Qgen     = 0.;
    buses[nameN].Pload    = 0.;
    buses[nameN].Qload    = 0.;
    buses[nameN].Gshunt   = 0.;
    buses[nameN].Bshunt   = 0.;
    buses[nameN].invX     = 0.;
    buses[nameN].bus_type = 1;
    buses[nameN].bus_name = id_name+".中性点"; //Matpower
    // ********************* Construction Complete ******************
    if (buses.count(name1) == 0 &&
        buses.count(name2) == 0 &&
        buses.count(name3) == 0)
      continue;
    
    ++count;
        if(Gmstr=="''"){
       Gm = 0.;
       //cout << "Gmstr" <<endl;
    }else{
       Gm = stod(Gmstr);
    }
    if(Bmstr=="\'\'"){
       Bm = 0.;
       //cout << "BBBmstr" <<endl;
    }else{
       Bm = stod(Bmstr); //string to double
    }
    if(abs(X1)<EPS) X1=0.001;
    if(abs(X2)<EPS) X2=0.001;
    if(abs(X3)<EPS) X3=0.001;
    // NEW ADD
    if(X1<=0) X1=-X1;
    if(X2<=0) X2=-X2;
    if(X3<=0) X3=-X3;  
    if(R1<=0) R1=-R1;
    if(R2<=0) R2=-R2;
    if(R3<=0) R3=-R3;

    /*if(I_off!=0){
    }*/   
    int MFILEPRINT_EN=0; //GraphDyn_util definition from EMSModel 
    int SC_SYS_INPUT=1; //GraphDyn_util definition from EMSModel
    if (MFILEPRINT_EN==0){
      if(nameN=="-5751"){
          cout << ">>> '-5751 neutral Joff='" << J_off << " J_island = " << J_island << endl;
        }  
      if(I_off==0 && I_island == 0){		
        if(SC_SYS_INPUT==1) tap1=-tap1;
        update_line(lines, nameN, name1, R1, X1, 0., Gm / 3., Bm / 3., -tap1, id_name+".高压端"); //dian ke yuan different tap sign //Matpower//
        update_line(lines, name1, nameN, R1, X1, 0., Gm / 3., Bm / 3., tap1, id_name+".高压端");
        if(nameN=="-5751"){
          cout << ">>> '-5751neutral R = '" << R1 << " X = " << X1 << endl;
        }
      }
      if(K_off==0 && K_island == 0){	
        if(SC_SYS_INPUT==1) tap2=-tap2;
        update_line(lines, nameN, name2, R2, X2, 0., Gm / 3., Bm / 3., -tap2, id_name+".中压端");//Matpower//
        update_line(lines, name2, nameN, R2, X2, 0., Gm / 3., Bm / 3., tap2, id_name+".中压端"); 
      }
      if(J_off==0 && J_island == 0){	
        if(SC_SYS_INPUT==1) tap3=-tap3;
        update_line(lines, nameN, name3, R3, X3, 0., Gm / 3., Bm / 3., -tap3, id_name+".低压端");//Matpower//
        update_line(lines, name3, nameN, R3, X3, 0., Gm / 3., Bm / 3., tap3, id_name+".低压端"); 
        if(nameN=="-5751")
        cout << ">>> '-5751 neutral R = '" << R3 << " X = " << X3 << endl;
      }
    }else {
      // if(I_off==0 && I_island == 0) update_line(lines, name1, nameN, R1, X1, 0., Gm / 3., Bm / 3., tap1, id_name+".高压端");
      // if(K_off==0 && K_island == 0) update_line(lines, name2, nameN, R2, X2, 0., Gm / 3., Bm / 3., tap2, id_name+".中压端");
      // if(J_off==0 && J_island == 0) update_line(lines, name3, nameN, R3, X3, 0., Gm / 3., Bm / 3., tap3, id_name+".低压端"); 
      if(I_off==0 && I_island == 0) update_line(lines, nameN, name1, R1, X1, 0., Gm / 3., Bm / 3., tap1, id_name+".高压端");//added on 6/23/2020 to 以前没有-号
      if(K_off==0 && K_island == 0) update_line(lines, nameN, name2, R2, X2, 0., Gm / 3., Bm / 3., tap2, id_name+".中压端");
      if(J_off==0 && J_island == 0) update_line(lines, nameN, name3, R3, X3, 0., Gm / 3., Bm / 3., tap3, id_name+".低压端"); 
    }
  }
  cout << "Three winding transformer count: " << count << endl;
}
//DO NOT NEED TO CHANGE RIGHT NOW//
void read_EPRI_GEN_data(unordered_map<int, EPRI_GEN_DATA> &all_gen, string path) {
  io::CSVReader<18> in(path);
  in.read_header(io::ignore_extra_column,
                 "Par_No", "XD", "XDP", "XDPP", "Xq",
                 "XqP", "XqPP", "X2", "Ra", "TD0P",
                 "TD0PP", "Tq0P", "Tq0PP", "TJ", "a",
                 "b", "n", "D_Sec"
                 );
  int Par_No;
  real__t Xd, Xdp, Xdpp, Xq, Xqp, Xqpp, X2, Ra, Td0p, Td0pp, Tq0p, Tq0pp, TJ, a, b, n, D;
  while(in.read_row(Par_No, Xd, Xdp, Xdpp, Xq,
                    Xqp, Xqpp, X2, Ra, Td0p,
                    Td0pp, Tq0p, Tq0pp, TJ, a,
                    b, n, D)){
    all_gen[Par_No].Xd    = Xd;
    all_gen[Par_No].Xdp   = Xdp < EPS  ? 0.0001 : Xdp;
    all_gen[Par_No].Xdpp  = Xdpp < EPS ? 0.0001 : Xdpp;
    all_gen[Par_No].Xq    = Xq;
    all_gen[Par_No].Xqp   = Xqp < EPS  ? 0.0001 : Xqp;
    all_gen[Par_No].Xqpp  = Xqpp < EPS ? 0.0001 : Xqpp;
    all_gen[Par_No].X2    = X2;
    all_gen[Par_No].Ra    = Ra;
    all_gen[Par_No].Td0p  = Td0p;
    all_gen[Par_No].Td0pp = Td0pp;
    all_gen[Par_No].Tq0p  = Tq0p;
    all_gen[Par_No].Tq0pp = Tq0pp;
    all_gen[Par_No].TJ    = TJ;
    all_gen[Par_No].a     = a;
    all_gen[Par_No].b     = b;
    all_gen[Par_No].n     = n;
    all_gen[Par_No].D     = D+500.;//2.;//
  }
}

void read_EPRI_GOV_I_data(unordered_map<int, EPRI_GOV_I_DATA> &all_gov_1, string path) {
  io::CSVReader<16> in(path);
  in.read_header(io::ignore_extra_column,
                 "Par_No", "Type", "Regu_Ratio", "T_Service", "Dead_Range",
                 "Valve_Max", "Valve_Min", "Gate_Max", "Gate_Min", "T0_Steam",
                 "Tw_Hydro", "K_Overheat", "T_Overheat", "K_Measure", "K_Sfeed",
                 "T_Sfeed"
                 );
  int Par_No, gov_type;
  real__t K_delta, TS, dead_band_tol, sigma_Max, sigma_Min, mu_Max, mu_Min, T0, TW, alpha, TRH, Ki, Kbeta, Ti;
  while(in.read_row(Par_No, gov_type, K_delta, TS, dead_band_tol,
                    sigma_Max, sigma_Min, mu_Max, mu_Min, T0,
                    TW, alpha, TRH, Ki, Kbeta,
                    Ti)){
    all_gov_1[Par_No].gov_type      = gov_type;
    all_gov_1[Par_No].K_delta       = K_delta;
    all_gov_1[Par_No].TS            = TS;
    all_gov_1[Par_No].dead_band_tol = dead_band_tol;
    all_gov_1[Par_No].sigma_Max     = sigma_Max;
    all_gov_1[Par_No].sigma_Min     = sigma_Min;
    all_gov_1[Par_No].mu_Max        = mu_Max;
    all_gov_1[Par_No].mu_Min        = mu_Min;
    all_gov_1[Par_No].T0            = T0;
    all_gov_1[Par_No].TW            = TW;
    all_gov_1[Par_No].alpha         = alpha;
    all_gov_1[Par_No].TRH           = TRH;
    all_gov_1[Par_No].Ki            = Ki;
    all_gov_1[Par_No].Kbeta         = Kbeta;
    all_gov_1[Par_No].Ti            = Ti;
  }
}

void read_EPRI_GOV_II_data(unordered_map<int, EPRI_GOV_II_DATA> &all_gov_2, string path) {
  io::CSVReader<18> in(path);
  in.read_header(io::ignore_extra_column,
                 "Par_No", "K_REV", "DEATH", "TR", "TB",
                 "Topen", "TC", "VELC", "VELO", "PMAX",
                 "PMIN", "TCH", "FHP", "TRH", "FIP",
                 "TCO", "FLP", "LND"
                 );
  int Par_No;
  real__t K, dead_band_tol, Tr, Tb, TO, TC, VEL_Close, VEL_Open, U_Max, U_Min, Tch, Fhp, Trh, Fip, Tco, Flp, lambda;
  while(in.read_row(Par_No, K, dead_band_tol, Tr, Tb,
                    TO, TC, VEL_Close, VEL_Open, U_Max,
                    U_Min, Tch, Fhp, Trh, Fip,
                    Tco, Flp, lambda)){
    all_gov_2[Par_No].K             = K;
    all_gov_2[Par_No].dead_band_tol = dead_band_tol;
    all_gov_2[Par_No].Tr            = Tr;
    all_gov_2[Par_No].Tb            = Tb;
    all_gov_2[Par_No].TO            = TO;
    all_gov_2[Par_No].TC            = TC;
    all_gov_2[Par_No].VEL_Close     = VEL_Close;
    all_gov_2[Par_No].VEL_Open      = VEL_Open;
    all_gov_2[Par_No].U_Max         = U_Max;
    all_gov_2[Par_No].U_Min         = U_Min;
    
    all_gov_2[Par_No].steam.Tch    = Tch;
    all_gov_2[Par_No].steam.Fhp    = Fhp;
    all_gov_2[Par_No].steam.Trh    = Trh;
    all_gov_2[Par_No].steam.Fip    = Fip;
    all_gov_2[Par_No].steam.Tco    = Tco;
    all_gov_2[Par_No].steam.Flp    = Flp;
    all_gov_2[Par_No].steam.lambda = lambda;
  }
}

void read_EPRI_GOV_III_data(unordered_map<int, EPRI_GOV_III_DATA> &all_gov_3, string path) {
  io::CSVReader<44> in(path);
  in.read_header(io::ignore_extra_column,
                 "Par_No",
                 "DEATH", "K_REV", "T1",
                 "KMLD", // control 1 //负荷自动开关
                 "KP1", "KD1", "KI1", "IMAX1", "IMIN1", "PMAX1", "PMIN1", // pid load
                 "KMLDF", //control 2
//                 "FFMAX", "FFMIN",  //一次调频负荷上下限, not used
                 "KMAP",  // control 3
                 "KP2", "KD2", "KI2", "IMAX2", "IMIN2", "PMAX2", "PMIN2", "CMAX", "CMIN", // pid pressure
                 "TC", "Topen", "VELC", "VELO", "PMAX", "PMIN", "TF",  //ehs part 1
                 "KP", "KD", "KI", "IMAX", "IMIN", "PIDMAX", "PIDMIN", //ehs part 2 -- PID
                 "TCH", "FHP", "TRH", "FIP", "TCO", "FLP", "LND"  // steam
//                 , "TSH", "TD", "TW", "K_MP", "TDLY", "VMAX", "VMIN", "T4" // not used
                 );
  int Par_No;
  real__t dead_band_tol, K1, T1;
  int load_pid_control;
  real__t Kp1, Kd1, Ki1, I_Max1, I_Min1, PID_Max1, PID_Min1;
  int load_feed_control, pressure_pid_control;
  real__t Kp2, Kd2, Ki2, I_Max2, I_Min2, PID_Max2, PID_Min2, CON_Max, CON_Min;
  real__t TC, TO, VEL_Close, VEL_Open, P_Max, P_Min, T2, Kp, Kd, Ki, I_Max, I_Min, PID_Max, PID_Min;
  real__t Tch, Fhp, Trh, Fip, Tco, Flp, lambda;
  
  while(in.read_row(Par_No, dead_band_tol, K1, T1, load_pid_control,
                    Kp1, Kd1, Ki1, I_Max1, I_Min1, PID_Max1, PID_Min1, load_feed_control, pressure_pid_control,
                    Kp2, Kd2, Ki2, I_Max2, I_Min2, PID_Max2, PID_Min2, CON_Max, CON_Min,
                    TC, TO, VEL_Close, VEL_Open, P_Max, P_Min, T2,
                    Kp, Kd, Ki, I_Max, I_Min, PID_Max, PID_Min,
                    Tch, Fhp, Trh, Fip, Tco, Flp, lambda)){
    all_gov_3[Par_No].dead_band_tol        = dead_band_tol;
    all_gov_3[Par_No].K1                   = K1;
    all_gov_3[Par_No].T1                   = T1;
    all_gov_3[Par_No].load_pid_control     = load_pid_control;
    all_gov_3[Par_No].load_feed_control    = load_feed_control;
    all_gov_3[Par_No].pressure_pid_control = pressure_pid_control;
    
    all_gov_3[Par_No].pid_load.Kp      = Kp1;
    all_gov_3[Par_No].pid_load.Kd      = Kd1;
    all_gov_3[Par_No].pid_load.Ki      = Ki1;
    all_gov_3[Par_No].pid_load.I_Max   = I_Max1;
    all_gov_3[Par_No].pid_load.I_Min   = I_Min1;
    all_gov_3[Par_No].pid_load.PID_Max = PID_Max1;
    all_gov_3[Par_No].pid_load.PID_Min = PID_Min1;
    
    all_gov_3[Par_No].pid_pressure.Kp      = Kp2;
    all_gov_3[Par_No].pid_pressure.Kd      = Kd2;
    all_gov_3[Par_No].pid_pressure.Ki      = Ki2;
    all_gov_3[Par_No].pid_pressure.I_Max   = I_Max2;
    all_gov_3[Par_No].pid_pressure.I_Min   = I_Min2;
    all_gov_3[Par_No].pid_pressure.PID_Max = PID_Max2;
    all_gov_3[Par_No].pid_pressure.PID_Min = PID_Min2;
    
    all_gov_3[Par_No].CON_Max = CON_Max;
    all_gov_3[Par_No].CON_Min = CON_Min;
    
    all_gov_3[Par_No].ehs.TC          = TC;
    all_gov_3[Par_No].ehs.TO          = TO;
    all_gov_3[Par_No].ehs.VEL_Open    = VEL_Open;
    all_gov_3[Par_No].ehs.VEL_Close   = VEL_Close;
    all_gov_3[Par_No].ehs.P_Max       = P_Max;
    all_gov_3[Par_No].ehs.P_Min       = P_Min;
    all_gov_3[Par_No].ehs.T2          = T2;
    all_gov_3[Par_No].ehs.pid.Kp      = Kp;
    all_gov_3[Par_No].ehs.pid.Kd      = Kd;
    all_gov_3[Par_No].ehs.pid.Ki      = Ki;
    all_gov_3[Par_No].ehs.pid.I_Max   = I_Max;
    all_gov_3[Par_No].ehs.pid.I_Min   = I_Min;
    all_gov_3[Par_No].ehs.pid.PID_Max = PID_Max;
    all_gov_3[Par_No].ehs.pid.PID_Min = PID_Min;
    
    all_gov_3[Par_No].steam.Tch    = Tch;
    all_gov_3[Par_No].steam.Fhp    = Fhp;
    all_gov_3[Par_No].steam.Trh    = Trh;
    all_gov_3[Par_No].steam.Fip    = Fip;
    all_gov_3[Par_No].steam.Tco    = Tco;
    all_gov_3[Par_No].steam.Flp    = Flp;
    all_gov_3[Par_No].steam.lambda = lambda;
  }
}

void read_EPRI_GOV_IV_data(unordered_map<int, EPRI_GOV_IV_DATA> &all_gov_4, string path) {
  io::CSVReader<34> in(path);
  in.read_header(io::ignore_extra_column,
                 "Par_No", "T1", "DEATH", "K_REV",
                 "KMCT",  // control
                 "KP1", "KD1", "KI1", "IMAX1", "IMIN1", "PMAX1", "PMIN1",  // pid
                 "K2",
                 //"FFMAX", "FFMIN",  // not used
                 "TC", "Topen", "VELC", "VELO", "PMAX", "PMIN", "TF",  // ehs part 1
                 "KP", "KD", "KI", "IMAX", "IMIN", "PIDMAX", "PIDMIN", // ehs part 2 -- PID
                 "TCH", "FHP", "TRH", "FIP", "TCO", "FLP", "LND"  // steam
                 //, "TSH", "TD", "TW", "K_MP", "TDLY", "VMAX", "VMIN", "T4" // not used
                 );
  int Par_No;
  real__t T1, dead_band_tol, K;
  int control;
  real__t Kp1, Kd1, Ki1, I_Max1, I_Min1, PID_Max1, PID_Min1;
  real__t K2;
  real__t TC, TO, VEL_Close, VEL_Open, P_Max, P_Min, T2, Kp, Kd, Ki, I_Max, I_Min, PID_Max, PID_Min;
  real__t Tch, Fhp, Trh, Fip, Tco, Flp, lambda;
  
  while(in.read_row(Par_No, T1, dead_band_tol, K, control,
                    Kp1, Kd1, Ki1, I_Max1, I_Min1, PID_Max1, PID_Min1, K2,
                    TC, TO, VEL_Close, VEL_Open, P_Max, P_Min, T2,
                    Kp, Kd, Ki, I_Max, I_Min, PID_Max, PID_Min,
                    Tch, Fhp, Trh, Fip, Tco, Flp, lambda)){
    all_gov_4[Par_No].T1             = T1;
    all_gov_4[Par_No].dead_band_tol  = dead_band_tol;
    all_gov_4[Par_No].K              = K;
    all_gov_4[Par_No].control_option = control;
    all_gov_4[Par_No].K2             = K2;

    all_gov_4[Par_No].pid.Kp      = Kp1;
    all_gov_4[Par_No].pid.Kd      = Kd1;
    all_gov_4[Par_No].pid.Ki      = Ki1;
    all_gov_4[Par_No].pid.I_Max   = I_Max1;
    all_gov_4[Par_No].pid.I_Min   = I_Min1;
    all_gov_4[Par_No].pid.PID_Max = PID_Max1;
    all_gov_4[Par_No].pid.PID_Min = PID_Min1;
    
    all_gov_4[Par_No].ehs.TC          = TC;
    all_gov_4[Par_No].ehs.TO          = TO;
    all_gov_4[Par_No].ehs.VEL_Open    = VEL_Open;
    all_gov_4[Par_No].ehs.VEL_Close   = VEL_Close;
    all_gov_4[Par_No].ehs.P_Max       = P_Max;
    all_gov_4[Par_No].ehs.P_Min       = P_Min;
    all_gov_4[Par_No].ehs.T2          = T2;
    all_gov_4[Par_No].ehs.pid.Kp      = Kp;
    all_gov_4[Par_No].ehs.pid.Kd      = Kd;
    all_gov_4[Par_No].ehs.pid.Ki      = Ki;
    all_gov_4[Par_No].ehs.pid.I_Max   = I_Max;
    all_gov_4[Par_No].ehs.pid.I_Min   = I_Min;
    all_gov_4[Par_No].ehs.pid.PID_Max = PID_Max;
    all_gov_4[Par_No].ehs.pid.PID_Min = PID_Min;

    all_gov_4[Par_No].steam.Tch    = Tch;
    all_gov_4[Par_No].steam.Fhp    = Fhp;
    all_gov_4[Par_No].steam.Trh    = Trh;
    all_gov_4[Par_No].steam.Fip    = Fip;
    all_gov_4[Par_No].steam.Tco    = Tco;
    all_gov_4[Par_No].steam.Flp    = Flp;
    all_gov_4[Par_No].steam.lambda = lambda;
  }
}

void read_EPRI_GOV_V_data(unordered_map<int, EPRI_GOV_V_DATA> &all_gov_5, string path) {
  io::CSVReader<35> in(path);
  in.read_header(io::ignore_extra_column,
                 "Par_No", "T_LB", "DEATH", "KMLB", "K1", "K2", "KMCT",
                 "KP1", "KD1", "KI1", "IMAX1", "IMIN1", "PMAX1", "PMIN1",  // pid
                 //"FFMAX", "FFMIN",  // not used
                 "TC", "Topen", "VELC", "VELO", "PMAX", "PMIN", "TF",   // ehs part 1
                 "KP", "KD", "KI", "IMAX", "IMIN", "PIDMAX", "PIDMIN",  // ehs part 2 -- PID
                 "TCH", "FHP", "TRH", "FIP", "TCO", "FLP", "LND" // steam
                 //, "TSH", "TD", "TW", "K_MP", "TDLY", "VMAX", "VMIN", "T4" // not used
                 );
  int Par_No;
  real__t T1, dead_band_tol;
  int control1;
  real__t K1, K2;
  int control2;
  real__t Kp1, Kd1, Ki1, I_Max1, I_Min1, PID_Max1, PID_Min1;
  real__t TC, TO, VEL_Close, VEL_Open, P_Max, P_Min, T2, Kp, Kd, Ki, I_Max, I_Min, PID_Max, PID_Min;
  real__t Tch, Fhp, Trh, Fip, Tco, Flp, lambda;
  
  while(in.read_row(Par_No, T1, dead_band_tol, control1, K1, K2, control2,
                    Kp1, Kd1, Ki1, I_Max1, I_Min1, PID_Max1, PID_Min1,
                    TC, TO, VEL_Close, VEL_Open, P_Max, P_Min, T2,
                    Kp, Kd, Ki, I_Max, I_Min, PID_Max, PID_Min,
                    Tch, Fhp, Trh, Fip, Tco, Flp, lambda)){
    all_gov_5[Par_No].T1             = T1;
    all_gov_5[Par_No].dead_band_tol  = dead_band_tol;
    all_gov_5[Par_No].control_speed  = control1;
    all_gov_5[Par_No].K1             = K1;
    all_gov_5[Par_No].K2             = K2;
    all_gov_5[Par_No].control_option = control2;
    
    all_gov_5[Par_No].pid.Kp      = Kp1;
    all_gov_5[Par_No].pid.Kd      = Kd1;
    all_gov_5[Par_No].pid.Ki      = Ki1;
    all_gov_5[Par_No].pid.I_Max   = I_Max1;
    all_gov_5[Par_No].pid.I_Min   = I_Min1;
    all_gov_5[Par_No].pid.PID_Max = PID_Max1;
    all_gov_5[Par_No].pid.PID_Min = PID_Min1;
    
    all_gov_5[Par_No].ehs.TC          = TC;
    all_gov_5[Par_No].ehs.TO          = TO;
    all_gov_5[Par_No].ehs.VEL_Open    = VEL_Open;
    all_gov_5[Par_No].ehs.VEL_Close   = VEL_Close;
    all_gov_5[Par_No].ehs.P_Max       = P_Max;
    all_gov_5[Par_No].ehs.P_Min       = P_Min;
    all_gov_5[Par_No].ehs.T2          = T2;
    all_gov_5[Par_No].ehs.pid.Kp      = Kp;
    all_gov_5[Par_No].ehs.pid.Kd      = Kd;
    all_gov_5[Par_No].ehs.pid.Ki      = Ki;
    all_gov_5[Par_No].ehs.pid.I_Max   = I_Max;
    all_gov_5[Par_No].ehs.pid.I_Min   = I_Min;
    all_gov_5[Par_No].ehs.pid.PID_Max = PID_Max;
    all_gov_5[Par_No].ehs.pid.PID_Min = PID_Min;
    
    all_gov_5[Par_No].steam.Tch    = Tch;
    all_gov_5[Par_No].steam.Fhp    = Fhp;
    all_gov_5[Par_No].steam.Trh    = Trh;
    all_gov_5[Par_No].steam.Fip    = Fip;
    all_gov_5[Par_No].steam.Tco    = Tco;
    all_gov_5[Par_No].steam.Flp    = Flp;
    all_gov_5[Par_No].steam.lambda = lambda;
  }
}

void read_EPRI_GOV_VII_data(unordered_map<int, EPRI_GOV_VII_DATA> &all_gov_7, string path) {
  io::CSVReader<53> in(path);
  in.read_header(io::ignore_extra_column,
                 "Par_No",
                 "TYPE1", "TYPE2", "TYPE3",
                 "T1", "TR1", "DB1", "DB1MAX", "DB1MIN", "KW",
                 "T2", "TR2", "EP", "DB2", "DB2MAX", "DB2MIN",
                 "T3", "TR3", "BP", "DB3", "DB3MAX", "DB3MIN",
                 "KP", "KD", "KI", "TD", "INTMAX", "INTMIN", "PIDMAX", "PIDMIN",  // PID
                 "YC",
                 "Tclose", "Topen", "VELC", "VELO", "PMAX", "PMIN", "TF",  // ehs 1
                 "KP1", "KD1", "KI1", "INTMAX1", "INTMIN1", "PIDMAX1", "PIDMIN1", "T4", // ehs 2
                 "Tw", "Kwtype",
                 // "At", "Tr", "YNL", "DBWP", "DBWN",       // hydro machine
                 "pchoice", "promax1", "promin1", /*"promax2", "promin2",*/ "ratelimp", "ratelimn" // bounds
                 );
  int Par_No;
  int auto_manual_mode, additional_control, Y_control_input_option, pchoice, hydro_type;
  real__t T1, TR1, dead_band1, dead_band_Max1, dead_band_Min1, Kw;
  real__t T2, TR2, Ep, dead_band2, dead_band_Max2, dead_band_Min2;
  real__t T3, TR3, Bp, dead_band3, dead_band_Max3, dead_band_Min3;
  real__t Kp, Kd, Ki, Td, I_Max, I_Min, PID_Max, PID_Min, YC;
  real__t TC, TO, VEL_Close, VEL_Open, P_Max, P_Min, Tf, Kp1, Kd1, Ki1, I_Max1, I_Min1, PID_Max1, PID_Min1, T4;
  real__t Tw, PRO_Max, PRO_Min, Ratelimp, Ratelimn;
  
  while(in.read_row(Par_No, auto_manual_mode, additional_control, Y_control_input_option,
                    T1, TR1, dead_band1, dead_band_Max1, dead_band_Min1, Kw,
                    T2, TR2, Ep, dead_band2, dead_band_Max2, dead_band_Min2,
                    T3, TR3, Bp, dead_band3, dead_band_Max3, dead_band_Min3,
                    Kp, Kd, Ki, Td, I_Max, I_Min, PID_Max, PID_Min, YC,
                    TC, TO, VEL_Close, VEL_Open, P_Max, P_Min, Tf,
                    Kp1, Kd1, Ki1, I_Max1, I_Min1, PID_Max1, PID_Min1, T4,
                    Tw, hydro_type, pchoice, PRO_Max, PRO_Min, Ratelimp, Ratelimn)){
    all_gov_7[Par_No].auto_manual_mode                  = auto_manual_mode;
    all_gov_7[Par_No].additional_control                = additional_control;
    all_gov_7[Par_No].Y_control_input_option            = Y_control_input_option;
    all_gov_7[Par_No].additional_control_input_position = pchoice;
    all_gov_7[Par_No].hydro_type                        = hydro_type;

    all_gov_7[Par_No].T1             = T1;
    all_gov_7[Par_No].TR1            = TR1;
    all_gov_7[Par_No].dead_band1     = dead_band1;
    all_gov_7[Par_No].dead_band_Max1 = dead_band_Max1;
    all_gov_7[Par_No].dead_band_Min1 = dead_band_Min1;
    all_gov_7[Par_No].Kw             = Kw;
    
    all_gov_7[Par_No].T2             = T2;
    all_gov_7[Par_No].TR2            = TR2;
    all_gov_7[Par_No].Ep             = Ep;
    all_gov_7[Par_No].dead_band2     = dead_band2;
    all_gov_7[Par_No].dead_band_Max2 = dead_band_Max2;
    all_gov_7[Par_No].dead_band_Min2 = dead_band_Min2;

    all_gov_7[Par_No].T3             = T3;
    all_gov_7[Par_No].TR3            = TR3;
    all_gov_7[Par_No].Bp             = Bp;
    all_gov_7[Par_No].dead_band3     = dead_band3;
    all_gov_7[Par_No].dead_band_Max3 = dead_band_Max3;
    all_gov_7[Par_No].dead_band_Min3 = dead_band_Min3;
        
    all_gov_7[Par_No].Kp      = Kp;
    all_gov_7[Par_No].Kd      = Kd;
    all_gov_7[Par_No].Ki      = Ki;
    all_gov_7[Par_No].Td      = Td;
    all_gov_7[Par_No].I_Max   = I_Max;
    all_gov_7[Par_No].I_Min   = I_Min;
    all_gov_7[Par_No].PID_Max = PID_Max;
    all_gov_7[Par_No].PID_Min = PID_Min;
    
    all_gov_7[Par_No].ehs.TC          = TC;
    all_gov_7[Par_No].ehs.TO          = TO;
    all_gov_7[Par_No].ehs.VEL_Open    = VEL_Open;
    all_gov_7[Par_No].ehs.VEL_Close   = VEL_Close;
    all_gov_7[Par_No].ehs.P_Max       = P_Max;
    all_gov_7[Par_No].ehs.P_Min       = P_Min;
    all_gov_7[Par_No].ehs.T2          = Tf;
    all_gov_7[Par_No].ehs.pid.Kp      = Kp1;
    all_gov_7[Par_No].ehs.pid.Kd      = Kd1;
    all_gov_7[Par_No].ehs.pid.Ki      = Ki1;
    all_gov_7[Par_No].ehs.pid.I_Max   = I_Max1;
    all_gov_7[Par_No].ehs.pid.I_Min   = I_Min1;
    all_gov_7[Par_No].ehs.pid.PID_Max = PID_Max1;
    all_gov_7[Par_No].ehs.pid.PID_Min = PID_Min1;
    all_gov_7[Par_No].T4              = T4;
    
    all_gov_7[Par_No].Tw       = Tw;
    all_gov_7[Par_No].PRO_Max  = PRO_Max;
    all_gov_7[Par_No].PRO_Min  = PRO_Min;
    all_gov_7[Par_No].Ratelimp = Ratelimp;
    all_gov_7[Par_No].Ratelimn = Ratelimn;
  }
}

/** 8型GOV仅实现开度模式 */
void read_EPRI_GOV_VIII_data(unordered_map<int, EPRI_GOV_VIII_DATA> &all_gov_8, string path) {
  io::CSVReader<30> in(path);
  in.read_header(io::ignore_extra_column,
                 "Par_No",
                 "T1", "TR1", "DB1P", "DB1N", "MAX1", "MIN1", "Rti1", "Rtd1",  //频率调节
                 "mode", //附加调节选择： 1-功率调节， 2-开度调节；
                 "NMC",  //bp反馈信号：1-开度，2-YPID；
                 "NPA",  //bp输入接入点：1-PID前，2-PID内积分前
                 "T3", "TR3", "BP", "DB3", "MAX3", "MIN3", //开度模式参数
                 "Kp2", "Ki2", "Kd2", "Td2", "INTMAX2", "INTMIN2", "PIDMAX2", "PIDMIN2", "KW", //开度模式PID
                 "Tw", "Rti0", "Rtd0"
                 );
  int Par_No;
  int additional_control, Y_control_input_option, additional_control_input_position;
  real__t T1, TR1, dead_band_p, dead_band_n, dead_band_Max1, dead_band_Min1, Kw, Rti0, Rtd0, Rti1, Rtd1;
  real__t T3, TR3, Bp, dead_band3, dead_band_Max3, dead_band_Min3;
  real__t Kp2, Ki2, Kd2, Td2, I_Max2, I_Min2, PID_Max2, PID_Min2;
  real__t Tw;

  while(in.read_row(Par_No, T1, TR1, dead_band_p, dead_band_n, dead_band_Max1, dead_band_Min1, Rti1, Rtd1,
                    additional_control, Y_control_input_option, additional_control_input_position,
                    T3, TR3, Bp, dead_band3, dead_band_Max3, dead_band_Min3,
                    Kp2, Ki2, Kd2, Td2, I_Max2, I_Min2, PID_Max2, PID_Min2, Kw,
                    Tw, Rti0, Rtd0)){
    all_gov_8[Par_No].additional_control                = additional_control;
    all_gov_8[Par_No].Y_control_input_option            = Y_control_input_option;
    all_gov_8[Par_No].additional_control_input_position = additional_control_input_position;
    
    all_gov_8[Par_No].T1             = T1;
    all_gov_8[Par_No].TR1            = TR1;
    all_gov_8[Par_No].dead_band_p    = dead_band_p;
    all_gov_8[Par_No].dead_band_n    = dead_band_n;
    all_gov_8[Par_No].dead_band_Max1 = dead_band_Max1;
    all_gov_8[Par_No].dead_band_Min1 = dead_band_Min1;
    all_gov_8[Par_No].Rti0           = Rti0;
    all_gov_8[Par_No].Rtd0           = Rtd0;
    all_gov_8[Par_No].Rti1           = Rti1;
    all_gov_8[Par_No].Rtd1           = Rtd1;
    
    all_gov_8[Par_No].T3             = T3;
    all_gov_8[Par_No].TR3            = TR3;
    all_gov_8[Par_No].Bp             = Bp;
    all_gov_8[Par_No].dead_band3     = dead_band3;
    all_gov_8[Par_No].dead_band_Max3 = dead_band_Max3;
    all_gov_8[Par_No].dead_band_Min3 = dead_band_Min3;

    all_gov_8[Par_No].Kp2      = Kp2;
    all_gov_8[Par_No].Kd2      = Kd2;
    all_gov_8[Par_No].Ki2      = Ki2;
    all_gov_8[Par_No].Td2      = Td2;
    all_gov_8[Par_No].I_Max2   = I_Max2;
    all_gov_8[Par_No].I_Min2   = I_Min2;
    all_gov_8[Par_No].PID_Max2 = PID_Max2;
    all_gov_8[Par_No].PID_Min2 = PID_Min2;

    all_gov_8[Par_No].Tw = Tw;
  }
}

void read_EPRI_GOV_IX_data(unordered_map<int, EPRI_GOV_IX_DATA> &all_gov_9, string path) {
  io::CSVReader<37> in(path);
  in.read_header(io::ignore_extra_column,
                 "Par_No", "T1", "DEATH", "K_REV",
                 "KMCT",  // control
                 "KP1", "KD1", "KI1", "IMAX1", "IMIN1", "PMAX1", "PMIN1",  // pid
                 "K2",
                 //"FFMAX", "FFMIN",  // not used
                 "TC", "Topen", "VELC", "VELO", "PMAX", "PMIN", "TF",  // ehs part 1
                 "KP", "KD", "KI", "IMAX", "IMIN", "PIDMAX", "PIDMIN", // ehs part 2 -- PID
                 "TCH", "FHP", "TRH", "FIP", "TCO", "FLP", "LND",  // steam
                 //, "TSH", "TD", "TW", "K_MP", "TDLY", "VMAX", "VMIN" // not used
                 "TDELAY1", "TDELAY2",  "TDELAY3"
                 );
  int Par_No;
  real__t T1, dead_band_tol, K;
  int control;
  real__t Kp1, Kd1, Ki1, I_Max1, I_Min1, PID_Max1, PID_Min1;
  real__t K2;
  real__t TC, TO, VEL_Close, VEL_Open, P_Max, P_Min, T2, Kp, Kd, Ki, I_Max, I_Min, PID_Max, PID_Min;
  real__t Tch, Fhp, Trh, Fip, Tco, Flp, lambda;
  real__t Tdelay1, Tdelay2, Tdelay3;
  
  while(in.read_row(Par_No, T1, dead_band_tol, K, control,
                    Kp1, Kd1, Ki1, I_Max1, I_Min1, PID_Max1, PID_Min1, K2,
                    TC, TO, VEL_Close, VEL_Open, P_Max, P_Min, T2,
                    Kp, Kd, Ki, I_Max, I_Min, PID_Max, PID_Min,
                    Tch, Fhp, Trh, Fip, Tco, Flp, lambda,
                    Tdelay1, Tdelay2, Tdelay3)){
    all_gov_9[Par_No].T1             = T1;
    all_gov_9[Par_No].dead_band_tol  = dead_band_tol;
    all_gov_9[Par_No].K              = K;
    all_gov_9[Par_No].control_option = control;
    all_gov_9[Par_No].K2             = K2;

    all_gov_9[Par_No].pid.Kp      = Kp1;
    all_gov_9[Par_No].pid.Kd      = Kd1;
    all_gov_9[Par_No].pid.Ki      = Ki1;
    all_gov_9[Par_No].pid.I_Max   = I_Max1;
    all_gov_9[Par_No].pid.I_Min   = I_Min1;
    all_gov_9[Par_No].pid.PID_Max = PID_Max1;
    all_gov_9[Par_No].pid.PID_Min = PID_Min1;
    
    all_gov_9[Par_No].ehs.TC          = TC;
    all_gov_9[Par_No].ehs.TO          = TO;
    all_gov_9[Par_No].ehs.VEL_Open    = VEL_Open;
    all_gov_9[Par_No].ehs.VEL_Close   = VEL_Close;
    all_gov_9[Par_No].ehs.P_Max       = P_Max;
    all_gov_9[Par_No].ehs.P_Min       = P_Min;
    all_gov_9[Par_No].ehs.T2          = T2;
    all_gov_9[Par_No].ehs.pid.Kp      = Kp;
    all_gov_9[Par_No].ehs.pid.Kd      = Kd;
    all_gov_9[Par_No].ehs.pid.Ki      = Ki;
    all_gov_9[Par_No].ehs.pid.I_Max   = I_Max;
    all_gov_9[Par_No].ehs.pid.I_Min   = I_Min;
    all_gov_9[Par_No].ehs.pid.PID_Max = PID_Max;
    all_gov_9[Par_No].ehs.pid.PID_Min = PID_Min;

    all_gov_9[Par_No].steam.Tch    = Tch;
    all_gov_9[Par_No].steam.Fhp    = Fhp;
    all_gov_9[Par_No].steam.Trh    = Trh;
    all_gov_9[Par_No].steam.Fip    = Fip;
    all_gov_9[Par_No].steam.Tco    = Tco;
    all_gov_9[Par_No].steam.Flp    = Flp;
    all_gov_9[Par_No].steam.lambda = lambda;
    
    all_gov_9[Par_No].Tdelay1 = Tdelay1;
    all_gov_9[Par_No].Tdelay2 = Tdelay2;
    all_gov_9[Par_No].Tdelay3 = Tdelay3;

  }
}

void read_EPRI_EXC_I_data(unordered_map<int, EPRI_EXC_I_DATA> &all_exc_1, string path) {
  io::CSVReader<10> in(path);
  in.read_header(io::ignore_extra_column,
                 "Par_No", "K_Measure", "T_Measure", "K_Amplify", "T_Amplify",
                 "K_Feed", "T_Feed", "T_Exciter", "Efd_Max", "Efd_Min"
                 );
  int Par_No;
  real__t Kr, Tr, Ka, Ta, Kf, Tf, Te, Efd_Max, Efd_Min;
  while(in.read_row(Par_No, Kr, Tr, Ka, Ta, Kf, Tf, Te, Efd_Max, Efd_Min)){
    all_exc_1[Par_No].Kr      = Kr;
    all_exc_1[Par_No].Tr      = Tr;
    all_exc_1[Par_No].Ka      = Ka;
    all_exc_1[Par_No].Ta      = Ta;
    all_exc_1[Par_No].Kf      = Kf;
    all_exc_1[Par_No].Tf      = Tf;
    all_exc_1[Par_No].Te      = Te;
    all_exc_1[Par_No].Efd_Max = Efd_Max;
    all_exc_1[Par_No].Efd_Min = Efd_Min;
  }
}

void read_EPRI_EXC_II_data(unordered_map<int, EPRI_EXC_II_DATA> &all_exc_2, string path) {
  io::CSVReader<17> in(path);
  in.read_header(io::ignore_extra_column,
                 "Par_No", "K_Measure", "T_Measure", "K_Type", "T1",
                 "T2", "T3", "T4", "K_Amplify", "T_Amplify",
                 "K_Voltage", "K_Current", "K_Arc", "Efd_Max", "Efd_Min",
                 "Vt_Max", "Vt_Min"
                 );
  int Par_No;
  real__t Kr, Tr, K2, T1, T2, T3, T4, Ka, Ta, Kpt, Kit, Ke, Efd_Max, Efd_Min, Vta, Vtb;
  while(in.read_row(Par_No, Kr, Tr, K2, T1,
                    T2, T3, T4, Ka, Ta,
                    Kpt, Kit, Ke, Efd_Max, Efd_Min,
                    Vta, Vtb)){
    all_exc_2[Par_No].Kr      = Kr;
    all_exc_2[Par_No].Tr      = Tr;
    all_exc_2[Par_No].K2      = K2;
    all_exc_2[Par_No].T1      = T1;
    all_exc_2[Par_No].T2      = T2;
    all_exc_2[Par_No].T3      = T3;
    all_exc_2[Par_No].T4      = T4;
    all_exc_2[Par_No].Ka      = Ka;
    all_exc_2[Par_No].Ta      = Ta;
    all_exc_2[Par_No].Kpt     = Kpt;
    all_exc_2[Par_No].Kit     = Kit;
    all_exc_2[Par_No].Ke      = Ke;
    all_exc_2[Par_No].Efd_Max = Efd_Max;
    all_exc_2[Par_No].Efd_Min = Efd_Min;
    all_exc_2[Par_No].Vta     = Vta;
    all_exc_2[Par_No].Vtb     = Vtb;
  }
}

void read_EPRI_EXC_III_TO_X_data(unordered_map<int, EPRI_EXC_III_TO_X_DATA> &all_exc_3_10, string path) {
  io::CSVReader<28> in(path);
  in.read_header(io::ignore_extra_column,
                 "Par_No", "xc_modul", "tr_measure", "k_dcplus", "kv_integ",
                 "t1_series", "t2_series", "t3_series", "t4_series", "ka_amplify",
                 "ta_amplify", "vamax", "vamin", "kf_paral", "tf_paral",
                 "kh1_feed", "kb_adjust", "t5_adjust", "vrmax", "vrmin",
                 "ke_exciter", "te_exciter", "vemax", "c1_exciter", "c2_exciter",
                 "kc_load", "kd_exciter", "efdmax"
                 //, "R_MODEL", "vmax1", "vmin1", "vmax2", "vmin2",
                 //"lelocation", "lemodel", "leparno", "oelocation", "oemodel",
                 //"oeparno", "oilocation", "oimodel", "oiparno"
                 );
  int Par_No;
  real__t Xc, Tr, K, Kv, T1, T2, T3, T4, Ka, Ta, Va_Max, Va_Min, Kf, Tf;
  real__t KH1, KB, T5, Vr_Max, Vr_Min, Ke, Te, Ve_Max, C1, C2, Kc, Kd, Efd_Max;
  while(in.read_row(Par_No, Xc, Tr, K, Kv,
                    T1, T2, T3, T4, Ka,
                    Ta, Va_Max, Va_Min, Kf, Tf,
                    KH1, KB, T5, Vr_Max, Vr_Min,
                    Ke, Te, Ve_Max, C1, C2,
                    Kc, Kd, Efd_Max)){
    all_exc_3_10[Par_No].Xc      = Xc;
    all_exc_3_10[Par_No].Tr      = Tr;
    all_exc_3_10[Par_No].K       = K;
    all_exc_3_10[Par_No].Kv      = Kv;
    all_exc_3_10[Par_No].T1      = T1;
    all_exc_3_10[Par_No].T2      = T2;
    all_exc_3_10[Par_No].T3      = T3;
    all_exc_3_10[Par_No].T4      = T4;
    all_exc_3_10[Par_No].Ka      = Ka;
    all_exc_3_10[Par_No].Ta      = Ta;
    all_exc_3_10[Par_No].Va_Max  = Va_Max;
    all_exc_3_10[Par_No].Va_Min  = Va_Min;
    all_exc_3_10[Par_No].Kf      = Kf;
    all_exc_3_10[Par_No].Tf      = Tf;
    all_exc_3_10[Par_No].KH1     = KH1;
    all_exc_3_10[Par_No].KB      = KB;
    all_exc_3_10[Par_No].T5      = T5;
    all_exc_3_10[Par_No].Vr_Max  = Vr_Max;
    all_exc_3_10[Par_No].Vr_Min  = Vr_Min;
    all_exc_3_10[Par_No].Ke      = Ke;
    all_exc_3_10[Par_No].Te      = Te;
    all_exc_3_10[Par_No].Ve_Max  = Ve_Max;
    all_exc_3_10[Par_No].C1      = C1;
    all_exc_3_10[Par_No].C2      = C2;
    all_exc_3_10[Par_No].Kc      = Kc;
    all_exc_3_10[Par_No].Kd      = Kd;
    all_exc_3_10[Par_No].Efd_Max = Efd_Max;
  }
}

void read_EPRI_EXC_XI_TO_XII_data(unordered_map<int, EPRI_EXC_XI_TO_XII_DATA> &all_exc_11_12, string path) {
  io::CSVReader<19> in(path);
  in.read_header(io::ignore_extra_column,
                 "Par_No", "xc_modul", "tr_measure", "k_dcplus", "kv_integ",
                 "t1_series", "t2_series", "t3_series", "t4_series", "ka_amplify",
                 "ta_amplify", "vamax", "vamin", "kf_paral", "tf_paral",
                 "vrmax", "vrmin", "kc_load", "Vs_POS"
                 //,"R_MODEL", "vmax1", "vmin1", "vmax2", "vmin2",
                 //"lelocation", "lemodel", "leparno", "oelocation", "oemodel",
                 //"oeparno", "oilocation", "oimodel", "oiparno"
                 );
  int Par_No;
  real__t Xc, Tr, K, Kv, T1, T2, T3, T4, Ka, Ta, Va_Max, Va_Min, Kf, Tf, Vr_Max, Vr_Min, Kc, Vs_Pos;
  while(in.read_row(Par_No, Xc, Tr, K, Kv,
                    T1, T2, T3, T4, Ka,
                    Ta, Va_Max, Va_Min, Kf, Tf,
                    Vr_Max, Vr_Min, Kc, Vs_Pos)){
    all_exc_11_12[Par_No].Xc     = Xc;
    all_exc_11_12[Par_No].Tr     = Tr;
    all_exc_11_12[Par_No].K      = K;
    all_exc_11_12[Par_No].Kv     = Kv;
    all_exc_11_12[Par_No].T1     = T1;
    all_exc_11_12[Par_No].T2     = T2;
    all_exc_11_12[Par_No].T3     = T3;
    all_exc_11_12[Par_No].T4     = T4;
    all_exc_11_12[Par_No].Ka     = Ka;
    all_exc_11_12[Par_No].Ta     = Ta;
    all_exc_11_12[Par_No].Va_Max = Va_Max;
    all_exc_11_12[Par_No].Va_Min = Va_Min;
    all_exc_11_12[Par_No].Kf     = Kf;
    all_exc_11_12[Par_No].Tf     = Tf;
    all_exc_11_12[Par_No].Vr_Max = Vr_Max;
    all_exc_11_12[Par_No].Vr_Min = Vr_Min;
    all_exc_11_12[Par_No].Kc     = Kc;
    all_exc_11_12[Par_No].Vs_Pos = Vs_Pos;
  }
}

void read_EPRI_PSS_I_data(unordered_map<int, EPRI_PSS_I_DATA> &all_pss_1, string path) {
  io::CSVReader<12> in(path);
  in.read_header(io::ignore_extra_column,
                 "Par_No", "K_Rotate", "K_Power", "K_Voltage", "K_Type",
                 "T_Washout", "T1_Shift", "T2_Shift", "T3_Shift", "T4_Shift",
                 "Vs_Max", "Vs_Min"
                 );
  int Par_No;
  real__t Kq1, Kq2, Kq3, K, Tq, T1e, T2e, T3e, T4e, VS_Max, VS_Min;
  while(in.read_row(Par_No, Kq1, Kq2, Kq3, K, Tq, T1e, T2e, T3e, T4e, VS_Max, VS_Min)){
    all_pss_1[Par_No].Kq1    = Kq1;
    all_pss_1[Par_No].Kq2    = Kq2;
    all_pss_1[Par_No].Kq3    = Kq3;
    all_pss_1[Par_No].K      = K;
    all_pss_1[Par_No].Tq     = Tq;
    all_pss_1[Par_No].T1e    = T1e;
    all_pss_1[Par_No].T2e    = T2e;
    all_pss_1[Par_No].T3e    = T3e;
    all_pss_1[Par_No].T4e    = T4e;
    all_pss_1[Par_No].VS_Max = VS_Max;
    all_pss_1[Par_No].VS_Min = VS_Min;
  }
}

void read_EPRI_PSS_II_data(unordered_map<int, EPRI_PSS_II_DATA> &all_pss_2, string path) {
  io::CSVReader<14> in(path);
  in.read_header(io::ignore_extra_column,
                 "Par_No", "k_rotate", "k_power", "k_voltage", //"t_measure",
                 "t_washout1", "t_washout2", "t1_shift", "t2_shift", "t3_shift",
                 "t4_shift", "t5_shift", "t6_shift", "vs_max", "vs_min"
                 );
  int Par_No;
  real__t Kw, Kp, Kt, TW1, TW2, T1, T2, T3, T4, T5, T6, VS_Max, VS_Min;
  while(in.read_row(Par_No, Kw, Kp, Kt, TW1, TW2, T1, T2, T3, T4, T5, T6, VS_Max, VS_Min)){
    all_pss_2[Par_No].Kw     = Kw;
    all_pss_2[Par_No].Kp     = Kp;
    all_pss_2[Par_No].Kt     = Kt;
    all_pss_2[Par_No].TW1    = TW1;
    all_pss_2[Par_No].TW2    = TW2;
    all_pss_2[Par_No].T1     = T1;
    all_pss_2[Par_No].T2     = T2;
    all_pss_2[Par_No].T3     = T3;
    all_pss_2[Par_No].T4     = T4;
    all_pss_2[Par_No].T5     = T5;
    all_pss_2[Par_No].T6     = T6;
    all_pss_2[Par_No].VS_Max = VS_Max;
    all_pss_2[Par_No].VS_Min = VS_Min;
  }
}

void read_EPRI_PSS_IV_VI_data(unordered_map<int, EPRI_PSS_IV_VI_DATA> &all_pss_4_6, string path) {
  io::CSVReader<24> in(path);
  in.read_header(io::ignore_extra_column,
                 "Par_No", "Trw", "T5", "T6", "T7",
                 "Kr", "Trp", "Tw", "Tw1", "Tw2",
                 "Ks", "T9", "T10", "T12", "Kw",
                 "Kp", "T1", "T2", "T13", "T14",
                 "T3", "T4", "VSMAX", "VSMIN"
                 );
  int Par_No;
  real__t Trw, T5, T6, T7, Kr, Trp, Tw, Tw1, Tw2, Ks, T9, T10, T12;
  real__t Kw, Kp, T1, T2, T13, T14, T3, T4, VS_Max, VS_Min;
  while(in.read_row(Par_No, Trw, T5, T6, T7, Kr, Trp, Tw, Tw1, Tw2, Ks, T9, T10, T12,
                    Kw, Kp, T1, T2, T13, T14, T3, T4, VS_Max, VS_Min)){
    all_pss_4_6[Par_No].Trw    = Trw;
    all_pss_4_6[Par_No].T5     = T5;
    all_pss_4_6[Par_No].T6     = T6;
    all_pss_4_6[Par_No].T7     = T7;
    all_pss_4_6[Par_No].Kr     = Kr;
    all_pss_4_6[Par_No].Trp    = Trp;
    all_pss_4_6[Par_No].Tw     = Tw;
    all_pss_4_6[Par_No].Tw1    = Tw1;
    all_pss_4_6[Par_No].Tw2    = Tw2;
    all_pss_4_6[Par_No].Ks     = Ks;
    all_pss_4_6[Par_No].T9     = T9;
    all_pss_4_6[Par_No].T10    = T10;
    all_pss_4_6[Par_No].T12    = T12;
    all_pss_4_6[Par_No].Kw     = Kw;
    all_pss_4_6[Par_No].Kp     = Kp;
    all_pss_4_6[Par_No].T1     = T1;
    all_pss_4_6[Par_No].T2     = T2;
    all_pss_4_6[Par_No].T13    = T13;
    all_pss_4_6[Par_No].T14    = T14;
    all_pss_4_6[Par_No].T3     = T3;
    all_pss_4_6[Par_No].T4     = T4;
    all_pss_4_6[Par_No].VS_Max = VS_Max;
    all_pss_4_6[Par_No].VS_Min = VS_Min;
  }
}

void read_EPRI_PSS_V_data(unordered_map<int, EPRI_PSS_V_DATA> &all_pss_5, string path) {
  io::CSVReader<14> in(path);
  in.read_header(io::ignore_extra_column,
                 "Par_No", "T1", "T2", "T3", "T4",
                 "T5", "A", "P", "K1", "T6",
                 "K2", "K", "VSMAX", "VSMIN"
                 );
  int Par_No;
  real__t T1, T2, T3, T4, T5, a, p, K1, T6, K2, K, VS_Max, VS_Min;
  while(in.read_row(Par_No, T1, T2, T3, T4, T5, a, p, K1, T6, K2, K, VS_Max, VS_Min)){
    all_pss_5[Par_No].T1     = T1;
    all_pss_5[Par_No].T2     = T2;
    all_pss_5[Par_No].T3     = T3;
    all_pss_5[Par_No].T4     = T4;
    all_pss_5[Par_No].T5     = T5;
    all_pss_5[Par_No].a      = a;
    all_pss_5[Par_No].p      = p;
    all_pss_5[Par_No].K1     = K1;
    all_pss_5[Par_No].T6     = T6;
    all_pss_5[Par_No].K2     = K2;
    all_pss_5[Par_No].K      = K;
    all_pss_5[Par_No].VS_Max = VS_Max;
    all_pss_5[Par_No].VS_Min = VS_Min;
  }
}

void read_EPRI_PSS_VIII_data(unordered_map<int, EPRI_PSS_VIII_DATA> &all_pss_8, string path) {
  io::CSVReader<11> in(path);
  in.read_header(io::ignore_extra_column,
                 "Par_No", "Kqv", "Tqv", "TQ1", "TQ1P",
                 "TQ2", "TQ2P", "TQ3", "TQ3P", "Vsmax",
                 "Vsmin"
                 );
  int Par_No;
  real__t Kqv, Tqv, Tq1, Tq1p, Tq2, Tq2p, Tq3, Tq3p, VS_Max, VS_Min;
  while(in.read_row(Par_No, Kqv, Tqv, Tq1, Tq1p, Tq2, Tq2p, Tq3, Tq3p, VS_Max, VS_Min)){
    all_pss_8[Par_No].Kqv    = Kqv;
    all_pss_8[Par_No].Tqv    = Tqv;
    all_pss_8[Par_No].Tq1    = Tq1;
    all_pss_8[Par_No].Tq1p   = Tq1p;
    all_pss_8[Par_No].Tq2    = Tq2;
    all_pss_8[Par_No].Tq2p   = Tq2p;
    all_pss_8[Par_No].Tq3    = Tq3;
    all_pss_8[Par_No].Tq3p   = Tq3p;
    all_pss_8[Par_No].VS_Max = VS_Max;
    all_pss_8[Par_No].VS_Min = VS_Min;
  }
}

}
