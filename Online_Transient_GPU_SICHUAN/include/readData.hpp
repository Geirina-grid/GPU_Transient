#ifndef readData_hpp
#define readData_hpp

namespace transient_analysis {

void read_bus_data(map<string, BUS> &buses, string path);
void read_fdpf_data(map<string, BUS> &buses, string path);
void read_load_data(map<string, BUS> &buses, string path);
void read_compensator_P_data(map<string, BUS> &buses, string path);
void read_generator_node_data(map<string, BUS> &buses, map<string, GENERATOR> &generators, string path);
void update_line(unordered_map<string, LINE> &lines, string from, string to, real__t R, real__t X,
                 real__t Bh, real__t Gm, real__t Bm, real__t tap);
void read_DC_line_data(map<string, BUS> &buses, unordered_map<string, LINE> &line, string path);
void read_AC_line_data(map<string, BUS> &buses, unordered_map<string, LINE> &line, string path);
void read_two_winding_transformer_data(map<string, BUS> &buses, unordered_map<string, LINE> &line, string path);
void read_three_winding_transformer_data(map<string, BUS> &buses, unordered_map<string, LINE> &line, string path);
void read_EPRI_GEN_data(unordered_map<int, EPRI_GEN_DATA> &all_gen, string path);
void read_EPRI_GOV_I_data(unordered_map<int, EPRI_GOV_I_DATA> &all_gov_1, string path);
void read_EPRI_GOV_II_data(unordered_map<int, EPRI_GOV_II_DATA> &all_gov_2, string path);
void read_EPRI_GOV_III_data(unordered_map<int, EPRI_GOV_III_DATA> &all_gov_3, string path);
void read_EPRI_GOV_IV_data(unordered_map<int, EPRI_GOV_IV_DATA> &all_gov_4, string path);
void read_EPRI_GOV_V_data(unordered_map<int, EPRI_GOV_V_DATA> &all_gov_5, string path);
void read_EPRI_GOV_VII_data(unordered_map<int, EPRI_GOV_VII_DATA> &all_gov_7, string path);
void read_EPRI_GOV_VIII_data(unordered_map<int, EPRI_GOV_VIII_DATA> &all_gov_8, string path);
void read_EPRI_GOV_IX_data(unordered_map<int, EPRI_GOV_IX_DATA> &all_gov_9, string path);
void read_EPRI_EXC_I_data(unordered_map<int, EPRI_EXC_I_DATA> &all_exc_1, string path);
void read_EPRI_EXC_II_data(unordered_map<int, EPRI_EXC_II_DATA> &all_exc_2, string path);
void read_EPRI_EXC_III_TO_X_data(unordered_map<int, EPRI_EXC_III_TO_X_DATA> &all_exc_3_10, string path);
void read_EPRI_EXC_XI_TO_XII_data(unordered_map<int, EPRI_EXC_XI_TO_XII_DATA> &all_exc_11_12, string path);
void read_EPRI_PSS_I_data(unordered_map<int, EPRI_PSS_I_DATA> &all_pss_1, string path);
void read_EPRI_PSS_II_data(unordered_map<int, EPRI_PSS_II_DATA> &all_pss_2, string path);
void read_EPRI_PSS_IV_VI_data(unordered_map<int, EPRI_PSS_IV_VI_DATA> &all_pss_4_6, string path);
void read_EPRI_PSS_V_data(unordered_map<int, EPRI_PSS_V_DATA> &all_pss_5, string path);
void read_EPRI_PSS_VIII_data(unordered_map<int, EPRI_PSS_VIII_DATA> &all_pss_8, string path);
}

#endif /* readData_h */
