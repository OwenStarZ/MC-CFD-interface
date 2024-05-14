import os
import shutil

# 获取当前脚本所在目录的绝对路径
script_dir = os.path.dirname(os.path.abspath(__file__))


# Part I. 读取 coupling.dat 中的网格信息
coupling_path = os.path.join(script_dir, 'coupling.dat')
with open(coupling_path, 'r') as coupling_file:
    lines = coupling_file.readlines()
    for i, line in enumerate(lines):
        if line.strip() == '#DIM':
            dim = lines[i+1].split()
        elif line.strip() == '#OUTDIMF':
            out_F_dim = lines[i+1].split()
        elif line.strip() == '#OUTDIMC':
            out_C_dim = lines[i+1].split()
        elif line.strip() == '#OUTDIMM':
            out_M_dim = lines[i+1].split()
        elif line.strip() == '#P_bndry in [cm] {xmin xmax ymin ymax zmin zmax}':
            P_bndry = lines[i+1].split()
        elif line.strip() == '#F_bndry in [m] {xmin xmax ymin ymax zmin zmax}':
            F_bndry = lines[i+1].split()
        elif line.strip() == '#C_bndry in [m] {xmin xmax ymin ymax zmin zmax}':
            C_bndry = lines[i+1].split()
        elif line.strip() == '#M_bndry in [m] {xmin xmax ymin ymax zmin zmax}':
            M_bndry = lines[i+1].split()



# Part II. 替换 udf.c 文件中的宏定义
udf_path = os.path.join(script_dir, 'libudf', 'src', 'udf.c')  # udf.c 文件的路径
temp_file_path = udf_path + '.tmp'  # 创建一个临时文件用于修改

with open(udf_path, 'r', encoding='utf-8') as udf_file, open(temp_file_path, 'w', encoding='utf-8') as temp_file:
    for line in udf_file:
        
        # 网格划分数定义
        
        if line.startswith('#define DIM0'):
            temp_file.write(f'#define DIM0 {dim[0]}\n')
        elif line.startswith('#define DIM1'):
            temp_file.write(f'#define DIM1 {dim[1]}\n')
        elif line.startswith('#define DIM2'):
            temp_file.write(f'#define DIM2 {dim[2]}\n')
        elif line.startswith('#define OUTDIMF0'):
            temp_file.write(f'#define OUTDIMF0 {out_F_dim[0]}\n')
        elif line.startswith('#define OUTDIMF1'):
            temp_file.write(f'#define OUTDIMF1 {out_F_dim[1]}\n')
        elif line.startswith('#define OUTDIMF2'):
            temp_file.write(f'#define OUTDIMF2 {out_F_dim[2]}\n')
        elif line.startswith('#define OUTDIMC0'):
            temp_file.write(f'#define OUTDIMC0 {out_C_dim[0]}\n')
        elif line.startswith('#define OUTDIMC1'):
            temp_file.write(f'#define OUTDIMC1 {out_C_dim[1]}\n')
        elif line.startswith('#define OUTDIMC2'):
            temp_file.write(f'#define OUTDIMC2 {out_C_dim[2]}\n')
        elif line.startswith('#define OUTDIMM0'):
            temp_file.write(f'#define OUTDIMM0 {out_M_dim[0]}\n')
        elif line.startswith('#define OUTDIMM1'):
            temp_file.write(f'#define OUTDIMM1 {out_M_dim[1]}\n')
        elif line.startswith('#define OUTDIMM2'):
            temp_file.write(f'#define OUTDIMM2 {out_M_dim[2]}\n')
            
        # 输出网格边界定义
        
        # 输出网格边界定义
        
        elif line.startswith('#define F_xmin'):
            temp_file.write(f'#define F_xmin {float(F_bndry[0])}\n')
        elif line.startswith('#define F_xmax'):
            temp_file.write(f'#define F_xmax {float(F_bndry[1])}\n')
        elif line.startswith('#define F_ymin'):
            temp_file.write(f'#define F_ymin {float(F_bndry[2])}\n')
        elif line.startswith('#define F_ymax'):
            temp_file.write(f'#define F_ymax {float(F_bndry[3])}\n')
        elif line.startswith('#define F_zmin'):
            temp_file.write(f'#define F_zmin {float(F_bndry[4])}\n')
        elif line.startswith('#define F_zmax'):
            temp_file.write(f'#define F_zmax {float(F_bndry[5])}\n')

        elif line.startswith('#define C_xmin'):
            temp_file.write(f'#define C_xmin {float(C_bndry[0])}\n')
        elif line.startswith('#define C_xmax'):
            temp_file.write(f'#define C_xmax {float(C_bndry[1])}\n')
        elif line.startswith('#define C_ymin'):
            temp_file.write(f'#define C_ymin {float(C_bndry[2])}\n')
        elif line.startswith('#define C_ymax'):
            temp_file.write(f'#define C_ymax {float(C_bndry[3])}\n')
        elif line.startswith('#define C_zmin'):
            temp_file.write(f'#define C_zmin {float(C_bndry[4])}\n')
        elif line.startswith('#define C_zmax'):
            temp_file.write(f'#define C_zmax {float(C_bndry[5])}\n')
            
        elif line.startswith('#define M_xmin'):
            temp_file.write(f'#define M_xmin {float(M_bndry[0])}\n')
        elif line.startswith('#define M_xmax'):
            temp_file.write(f'#define M_xmax {float(M_bndry[1])}\n')
        elif line.startswith('#define M_ymin'):
            temp_file.write(f'#define M_ymin {float(M_bndry[2])}\n')
        elif line.startswith('#define M_ymax'):
            temp_file.write(f'#define M_ymax {float(M_bndry[3])}\n')
        elif line.startswith('#define M_zmin'):
            temp_file.write(f'#define M_zmin {float(M_bndry[4])}\n')
        elif line.startswith('#define M_zmax'):
            temp_file.write(f'#define M_zmax {float(M_bndry[5])}\n')
            
        else:
            temp_file.write(line)

shutil.copyfile(temp_file_path, udf_path)  # 强制替换原始文件

os.remove(temp_file_path)  # 删除临时文件


# Part III. 替换 h5rw.cpp 文件中的宏定义
h5rw_path = os.path.join(script_dir, 'lib_h5rw', 'src', 'h5rw.cpp')  # h5rw.cpp 文件的路径
temp_file_path = h5rw_path + '.tmp'  # 创建一个临时文件用于修改

with open(h5rw_path, 'r', encoding='utf-8') as h5rw_file, open(temp_file_path, 'w', encoding='utf-8') as temp_file:
    for line in h5rw_file:
        
        # 网格划分数定义
        
        if line.startswith('#define DIM0'):
            temp_file.write(f'#define DIM0 {dim[0]}\n')
        elif line.startswith('#define DIM1'):
            temp_file.write(f'#define DIM1 {dim[1]}\n')
        elif line.startswith('#define DIM2'):
            temp_file.write(f'#define DIM2 {dim[2]}\n')
        elif line.startswith('#define OUTDIMF0'):
            temp_file.write(f'#define OUTDIMF0 {out_F_dim[0]}\n')
        elif line.startswith('#define OUTDIMF1'):
            temp_file.write(f'#define OUTDIMF1 {out_F_dim[1]}\n')
        elif line.startswith('#define OUTDIMF2'):
            temp_file.write(f'#define OUTDIMF2 {out_F_dim[2]}\n')
        elif line.startswith('#define OUTDIMC0'):
            temp_file.write(f'#define OUTDIMC0 {out_C_dim[0]}\n')
        elif line.startswith('#define OUTDIMC1'):
            temp_file.write(f'#define OUTDIMC1 {out_C_dim[1]}\n')
        elif line.startswith('#define OUTDIMC2'):
            temp_file.write(f'#define OUTDIMC2 {out_C_dim[2]}\n')
        elif line.startswith('#define OUTDIMM0'):
            temp_file.write(f'#define OUTDIMM0 {out_M_dim[0]}\n')
        elif line.startswith('#define OUTDIMM1'):
            temp_file.write(f'#define OUTDIMM1 {out_M_dim[1]}\n')
        elif line.startswith('#define OUTDIMM2'):
            temp_file.write(f'#define OUTDIMM2 {out_M_dim[2]}\n')
        else:
            temp_file.write(line)

shutil.copyfile(temp_file_path, h5rw_path)  # 强制替换原始文件

os.remove(temp_file_path)  # 删除临时文件


# Part IV. 替换 init.cpp 文件中的宏定义
init_path = os.path.join(script_dir, 'init.cpp')  # init.cpp 文件的路径
temp_file_path = init_path + '.tmp'  # 创建一个临时文件用于修改

with open(init_path, 'r', encoding='utf-8') as init_file, open(temp_file_path, 'w', encoding='utf-8') as temp_file:
    for line in init_file:

        # 网格划分数定义
        
        if line.startswith('#define OUTDIMF0'):
            temp_file.write(f'#define OUTDIMF0 {out_F_dim[0]}\n')
        elif line.startswith('#define OUTDIMF1'):
            temp_file.write(f'#define OUTDIMF1 {out_F_dim[1]}\n')
        elif line.startswith('#define OUTDIMF2'):
            temp_file.write(f'#define OUTDIMF2 {out_F_dim[2]}\n')
        elif line.startswith('#define OUTDIMC0'):
            temp_file.write(f'#define OUTDIMC0 {out_C_dim[0]}\n')
        elif line.startswith('#define OUTDIMC1'):
            temp_file.write(f'#define OUTDIMC1 {out_C_dim[1]}\n')
        elif line.startswith('#define OUTDIMC2'):
            temp_file.write(f'#define OUTDIMC2 {out_C_dim[2]}\n')
        elif line.startswith('#define OUTDIMM0'):
            temp_file.write(f'#define OUTDIMM0 {out_M_dim[0]}\n')
        elif line.startswith('#define OUTDIMM1'):
            temp_file.write(f'#define OUTDIMM1 {out_M_dim[1]}\n')
        elif line.startswith('#define OUTDIMM2'):
            temp_file.write(f'#define OUTDIMM2 {out_M_dim[2]}\n')
            
        # 输出网格边界定义
        
        elif line.startswith('#define F_xmin'):
            temp_file.write(f'#define F_xmin {float(F_bndry[0])*100}\n')
        elif line.startswith('#define F_xmax'):
            temp_file.write(f'#define F_xmax {float(F_bndry[1])*100}\n')
        elif line.startswith('#define F_ymin'):
            temp_file.write(f'#define F_ymin {float(F_bndry[2])*100}\n')
        elif line.startswith('#define F_ymax'):
            temp_file.write(f'#define F_ymax {float(F_bndry[3])*100}\n')
        elif line.startswith('#define F_zmin'):
            temp_file.write(f'#define F_zmin {float(F_bndry[4])*100}\n')
        elif line.startswith('#define F_zmax'):
            temp_file.write(f'#define F_zmax {float(F_bndry[5])*100}\n')

        elif line.startswith('#define C_xmin'):
            temp_file.write(f'#define C_xmin {float(C_bndry[0])*100}\n')
        elif line.startswith('#define C_xmax'):
            temp_file.write(f'#define C_xmax {float(C_bndry[1])*100}\n')
        elif line.startswith('#define C_ymin'):
            temp_file.write(f'#define C_ymin {float(C_bndry[2])*100}\n')
        elif line.startswith('#define C_ymax'):
            temp_file.write(f'#define C_ymax {float(C_bndry[3])*100}\n')
        elif line.startswith('#define C_zmin'):
            temp_file.write(f'#define C_zmin {float(C_bndry[4])*100}\n')
        elif line.startswith('#define C_zmax'):
            temp_file.write(f'#define C_zmax {float(C_bndry[5])*100}\n')
            
        elif line.startswith('#define M_xmin'):
            temp_file.write(f'#define M_xmin {float(M_bndry[0])*100}\n')
        elif line.startswith('#define M_xmax'):
            temp_file.write(f'#define M_xmax {float(M_bndry[1])*100}\n')
        elif line.startswith('#define M_ymin'):
            temp_file.write(f'#define M_ymin {float(M_bndry[2])*100}\n')
        elif line.startswith('#define M_ymax'):
            temp_file.write(f'#define M_ymax {float(M_bndry[3])*100}\n')
        elif line.startswith('#define M_zmin'):
            temp_file.write(f'#define M_zmin {float(M_bndry[4])*100}\n')
        elif line.startswith('#define M_zmax'):
            temp_file.write(f'#define M_zmax {float(M_bndry[5])*100}\n')
            
        elif line.startswith('#define INIT_temp_fuel'):
            temp_file.write(f'#define INIT_temp_fuel {float(lines[6])}\n')
        elif line.startswith('#define INIT_temp_moderator'):
            temp_file.write(f'#define INIT_temp_moderator {float(lines[7])}\n')
        elif line.startswith('#define INIT_temp_coolant'):
            temp_file.write(f'#define INIT_temp_coolant {float(lines[8])}\n')
        elif line.startswith('#define INIT_r_coolant'):
            temp_file.write(f'#define INIT_r_coolant {float(lines[9])}\n')
        else:
            temp_file.write(line)

shutil.copyfile(temp_file_path, init_path)  # 强制替换原始文件

os.remove(temp_file_path)  # 删除临时文件


# Part V. 替换 .rmc 输入卡中的网格计数器定义
rmc_files = [filename for filename in os.listdir(script_dir) if filename.endswith('.rmc')]

if len(rmc_files) == 1:
    rmc_path = os.path.join(script_dir, rmc_files[0])  # .rmc 文件的路径
    print("The only existing .rmc file found:", rmc_path)
else:
    print("No or multiple .rmc files found in the directory.")

temp_file_path = rmc_path + '.tmp'  # 创建一个临时文件用于修改

with open(rmc_path, 'r', encoding='utf-8') as rmc_file, open(temp_file_path, 'w', encoding='utf-8') as temp_file:
    for line in rmc_file:
        
        # 网格划分数定义
        
        if line.startswith('                    Bound ='):
            temp_file.write(f'                    Bound = {P_bndry[0]} {P_bndry[1]} {P_bndry[2]} {P_bndry[3]} {P_bndry[4]} {P_bndry[5]}\n')
        elif line.startswith('                    Scope ='):
            temp_file.write(f'                    Scope = {dim[0]} {dim[1]} {dim[2]}\n')
        else:
            temp_file.write(line)

shutil.copyfile(temp_file_path, rmc_path)  # 强制替换原始文件

os.remove(temp_file_path)  # 删除临时文件